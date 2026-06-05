# Snakefile
import os
import json
import urllib.request
import urllib.parse
import csv



# Load configuration file variables
configfile: "config.yaml"

BASE_DIR = config["base_dir"]
PP_QC_DIR = config["pp_QC_dir"]
SRP_ID = config["srp_id"]
RUN_TABLE = os.path.join(BASE_DIR, "SraRunTable.csv")
#SRR_LIST_PATH = os.path.join(BASE_DIR, "SRR_Acc_List.txt")

os.makedirs(BASE_DIR, exist_ok=True)
os.makedirs(PP_QC_DIR, exist_ok=True)

if not os.path.exists(RUN_TABLE):
    print(f"SraRunTable.csv not found. Fetching from NCBI for {SRP_ID}...")
    
    search_params = {"db": "sra", "term": SRP_ID, "retmode": "json", "retmax": 500}
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(search_params)
    
    try:
        # Search SRA
        with urllib.request.urlopen(search_url) as response:
            id_list = json.loads(response.read().decode()).get("esearchresult", {}).get("idlist", [])
            
        if id_list:
            # Fetch Metadata
            fetch_params = {"db": "sra", "id": ",".join(id_list), "rettype": "runinfo", "retmode": "text"}
            fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(fetch_params)
            urllib.request.urlretrieve(fetch_url, RUN_TABLE)
            print("SraRunTable.csv downloaded successfully.")
    except Exception as e:
        print(f"Warning: Failed to retrieve SRA metadata: {e}")


# 1. Parse SraRunTable to dynamically establish SRR (Run) -> GSM (GEO Sample) mapping
# Dynamically map GSM -> list of R1 SRRs and GSM -> list of R2 SRRs based on spot length
gsm_to_download = {}
srrs_to_download = []

if os.path.exists(RUN_TABLE):
    with open(RUN_TABLE, 'r') as f:
        first_line = f.readline()
        delimiter = '\t' if '\t' in first_line else ','
    with open(RUN_TABLE, 'r') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        headers = reader.fieldnames
        
        run_col = next((h for h in headers if h.lower() in ['run', 'run_s']), None)
        gsm_col = next((h for h in headers if h.lower() in ['samplename']), None)
        
        
        if run_col and gsm_col:
            for row in reader:
                srr = row[run_col].strip()
                gsm = row[gsm_col].strip()
                
                if srr and gsm:
                     if gsm not in gsm_to_download:
                        gsm_to_download[gsm] = []
                     gsm_to_download[gsm].append(srr)
                     srrs_to_download.append(srr)
                    # Note: spot_len == 8 (Index reads) are ignored entirely

GSMS = sorted(list(set(gsm_to_download.keys())))
SRRS = sorted(list(set(srrs_to_download)))

# Target Rule: Defines final required outputs
rule all:
    input:
        os.path.join(BASE_DIR, "GSM_Acc_List.txt"),
        os.path.join(BASE_DIR, f"{SRP_ID}_processed.h5ad")


# Output is set to temp() so raw SRA fastqs are deleted once merged to GSMs, conserving storage.
# STEP 1: Download and combine SRRs per GSM, cleaning up raw files on-the-fly
# Output is set to temp() so the large combined fastqs are automatically deleted 
# as soon as STAR alignment completes.
rule download_and_concatenate_gsm:
    output:
        r1 = temp(os.path.join(BASE_DIR, "combined_fastq/{gsm}_1.fastq.gz")),
        r2 = temp(os.path.join(BASE_DIR, "combined_fastq/{gsm}_2.fastq.gz"))
    threads: 4
    conda:
        "scanpy_env"  # Activates your environment natively for the whole rule
    resources:
        mem_mb = 20000,
        runtime = 1440,
        download_slots = 1  # Strictly throttles downloads to 1 GSM at a time
    params:
        # FIX: Join with space to create a valid Bash sequence
        srrs_f = lambda wildcards: " ".join(gsm_to_download[wildcards.gsm]),
        temp_dir = os.path.join(BASE_DIR, "fastq/{gsm}_temp")
    shell:
        """
        # Ensure temporary and output folders exist
        mkdir -p {params.temp_dir}
        mkdir -p $(dirname {output.r1})

        # Truncate target files to ensure they are empty before appending
        > {output.r1}
        > {output.r2}

        # Download and append paired-end runs sequentially
        for srr in {params.srrs_f}; do
            # 1. Download SRA archive
            prefetch --max-size 100G -O {params.temp_dir} $srr
            
            # 2. Extract paired fastq files (_1.fastq and _2.fastq)
            fasterq-dump -e {threads} --outdir {params.temp_dir} {params.temp_dir}/$srr
            
            # 3. Compress both reads separately
            pigz -f -p {threads} {params.temp_dir}/${srr}_1.fastq
            pigz -f -p {threads} {params.temp_dir}/${srr}_2.fastq
            
            # 4. Append to target output (gzip format natively supports concatenation!)
            cat {params.temp_dir}/${srr}_1.fastq.gz >> {output.r1}
            cat {params.temp_dir}/${srr}_2.fastq.gz >> {output.r2}
            
            # 5. Clean up temporary files for this run immediately to conserve space
            rm -rf {params.temp_dir}/$srr {params.temp_dir}/${srr}_1.fastq.gz {params.temp_dir}/${srr}_2.fastq.gz
        done

        # Final cleanup of empty temp directory
        rm -rf {params.temp_dir}
        """




# STEP 2: Write unique GSMs list to a text file (Replaces combined_fastq.py final task)
rule generate_gsm_list:
    input:
        run_table = RUN_TABLE
    output:
        gsm_list = os.path.join(BASE_DIR, "GSM_Acc_List.txt")
    run:
        with open(output.gsm_list, 'w') as out_f:
            for gsm in GSMS:
                out_f.write(f"{gsm}\n")


# STEP 3: Align with STAR Solo (Keeps your exact read file arrangement and cell filters)
rule star_align:
    input:
        r1 = os.path.join(BASE_DIR, "combined_fastq/{gsm}_1.fastq.gz"),
        r2 = os.path.join(BASE_DIR, "combined_fastq/{gsm}_2.fastq.gz")
    output:
        bam = os.path.join(BASE_DIR, f"{SRP_ID}_STAR/{{gsm}}_Aligned.sortedByCoord.out.bam"),
        barcodes = os.path.join(BASE_DIR, f"{SRP_ID}_STAR/{{gsm}}_Solo.out/Gene/filtered/barcodes.tsv")
    threads: 12
    resources:
        mem_mb = 64000,
        runtime = 2400,
        star_slots = 1 
    params:
        genome_dir = config["genome_dir"],
        whitelist = config["solo_whitelist_v2"],
        prefix = os.path.join(BASE_DIR, f"{SRP_ID}_STAR/{{gsm}}_")
    shell:
        """
        # STAR alignment with R2 passed first and R1 passed second to match your config
        conda run -n scanpy_env STAR --runThreadN {threads} \
             --genomeDir {params.genome_dir} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.prefix} \
             --readFilesIn {input.r2} {input.r1} \
             --soloType CB_UMI_Simple \
             --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
             --soloBarcodeReadLength 0 \
             --soloCBwhitelist "{params.whitelist}" \
             --soloCellFilter EmptyDrops_CR \
             --soloStrand Forward \
             --outSAMattributes NH HI AS nM CB UB \
             --outSAMtype BAM SortedByCoordinate
        """


# STEP 4: Run Velocyto (Specifies high-RAM config mirroring your sbatch allocation)
rule velocyto:
    input:
        bam = os.path.join(BASE_DIR, f"{SRP_ID}_STAR/{{gsm}}_Aligned.sortedByCoord.out.bam"),
        barcodes = os.path.join(BASE_DIR, f"{SRP_ID}_STAR/{{gsm}}_Solo.out/Gene/filtered/barcodes.tsv")
    output:
        loom = os.path.join(BASE_DIR, f"{SRP_ID}_velocyto/{{gsm}}_velocyto/{{gsm}}_Aligned.sortedByCoord.out.loom")
    threads: 4
    resources:
        mem_mb = 200000, 
        runtime = 2880,
        velocyto_slots = 1  # <--- Assign 1 slot per Velocyto job
    params:
        out_dir = os.path.join(BASE_DIR, f"{SRP_ID}_velocyto/{{gsm}}_velocyto"),
        mask_file = config["mask_file"],
        gtf_file = config["gtf_file"]
    shell:
        """
        conda run -n scanpy_env velocyto run -b {input.barcodes} \
                     -o {params.out_dir} \
                     -m {params.mask_file} \
                     -e {wildcards.gsm}_Aligned.sortedByCoord.out \
                     -@ {threads} \
                     --samtools-memory 40000 \
                     {input.bam} \
                     {params.gtf_file}
        """


# STEP 5: Aggregation & QC Step (Combines all resolved loom outputs)
rule pp_qc:
    input:
        looms = expand(os.path.join(BASE_DIR, f"{SRP_ID}_velocyto/{{gsm}}_velocyto/{{gsm}}_Aligned.sortedByCoord.out.loom"), gsm=GSMS),
        script = f"data_pp_smk.py"
    output:
        h5ad = os.path.join(BASE_DIR, f"{SRP_ID}_processed.h5ad")
    threads: 4
    resources:
        mem_mb = 100000,
        runtime = 2400
    params:
        srp_id = SRP_ID,
        pp_qc_dir = PP_QC_DIR
    shell:
        """
        mv {input.script} ./{params.pp_qc_dir}/

        # Explicitly preload your environment's modern C++ library to resolve CXXABI errors
        export LD_PRELOAD=/home/hlinh/.conda/envs/scanpy_env/lib/libstdc++.so.6

        # Execute your script
        conda run -n scanpy_env python {input.script} --inputs {input.looms} --output {output.h5ad} --srp {params.srp_id}
        """