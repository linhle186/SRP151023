#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=208GB
#SBATCH --mail-user=gyu@genzentrum.lmu.de
#SBATCH --mail-type=fail
#SBATCH -o %j.log
#SBATCH -e %j.err
#SBATCH -J RNAvelo.pp.pipeline
#SBATCH -t 5-00:00:00

set -euo pipefail

# ===========================================================================
# 0. CONFIG  (mirrors config.yaml / config_prep.yaml)
# ===========================================================================
SRP_ID="SRP151023"
BASE_DIR="./data"
PP_QC_DIR="./pp_QC"

REF_DIR="references"
REF_HOMO="${REF_DIR}/homo"
REF_MUS="${REF_DIR}/mus"

# Active references for the main run = human + 10x v2 (matches config.yaml)
GENOME_DIR="${REF_HOMO}/star_index_homo"
GTF_FILE="${REF_HOMO}/Homo_sapiens.GRCh38.115.gtf"
MASK_FILE="${REF_HOMO}/hg38_repeat_masker.gtf"
WHITELIST_V2="${REF_DIR}/737K-august-2016.txt"
WHITELIST_V3="${REF_DIR}/3M-february-2018.txt"

# Thread counts kept identical to the Snakefile rules so the RAM math holds
# (velocyto: VELO_THREADS * samtools-memory must stay under the node's RAM).
STAR_THREADS="${STAR_THREADS:-12}"
VELO_THREADS="${VELO_THREADS:-4}"
DL_THREADS="${DL_THREADS:-4}"

# Snakefile_prep builds a mouse index too, but the human main run never uses it.
# Default: skip it (saves an index build). Set BUILD_MOUSE=true to build it.
BUILD_MOUSE="${BUILD_MOUSE:-false}"

echo "=========================================================="
echo "   SRP151023 PIPELINE"
echo "   $(date)   cwd=$(pwd)"
echo "=========================================================="

# ===========================================================================
# 1. CONDA ENV  (create from environment.yml if missing, then activate)
# ===========================================================================
source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda env list | grep -q "scanpy_env"; then
    echo "scanpy_env not found -> creating from environment.yml ..."
    conda env create -f environment.yml -n scanpy_env
fi
echo "Activating conda environment: scanpy_env"
conda activate scanpy_env

mkdir -p "$BASE_DIR" "$PP_QC_DIR"

# ===========================================================================
# 2. STAGE 1 -- REFERENCES (idempotent: each helper skips if output exists)
# ===========================================================================
dl_gunzip () {          # $1=url  $2=output (uncompressed)
    local url="$1" out="$2"
    if [ -f "$out" ]; then echo "  [skip] $out"; return; fi
    mkdir -p "$(dirname "$out")"
    wget -q -O "${out}.gz" "$url"
    gunzip -f "${out}.gz"
}

echo "=========================================================="
make_mask () {          # $1=url(rmsk.txt.gz)  $2=output.gtf
    local url="$1" out="$2"
    if [ -f "$out" ]; then echo "  [skip] $out"; return; fi
    mkdir -p "$(dirname "$out")"
    wget -q -O "${out}.rmsk.txt.gz" "$url"
    zcat "${out}.rmsk.txt.gz" \
      | awk 'BEGIN{OFS="\t"}{print $6,"RepeatMasker","exon",$7+1,$8,".",$10,".","gene_id \""$12"\"; transcript_id \""$12"\";"}' \
      > "$out"
    rm -f "${out}.rmsk.txt.gz"
}

build_index () {        # $1=fasta  $2=gtf  $3=index_dir
    local fasta="$1" gtf="$2" idx="$3"
    if [ -f "${idx}/Genome" ]; then echo "  [skip] STAR index $idx"; return; fi
    mkdir -p "$idx"
    STAR --runMode genomeGenerate \
         --runThreadN "$STAR_THREADS" \
         --genomeDir "$idx" \
         --genomeFastaFiles "$fasta" \
         --sjdbGTFfile "$gtf" \
         --sjdbOverhang 99
}

echo "----- STAGE 1: references -----"
dl_gunzip "http://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
          "${REF_HOMO}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
dl_gunzip "http://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz" "$GTF_FILE"
make_mask "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz" "$MASK_FILE"

# 10x barcode whitelists
if [ ! -f "$WHITELIST_V2" ]; then mkdir -p "$REF_DIR"; wget -q -O "$WHITELIST_V2" \
  "https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/737K-august-2016.txt"; fi
if [ ! -f "$WHITELIST_V3" ]; then wget -q -O "$WHITELIST_V3" \
  "https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/3M-february-2018.txt"; fi

build_index "${REF_HOMO}/Homo_sapiens.GRCh38.dna.primary_assembly.fa" "$GTF_FILE" "$GENOME_DIR"

if [ "$BUILD_MOUSE" = "true" ]; then
    dl_gunzip "http://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" \
              "${REF_MUS}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
    dl_gunzip "http://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz" \
              "${REF_MUS}/Mus_musculus.GRCm39.115.gtf"
    make_mask "http://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/rmsk.txt.gz" "${REF_MUS}/mm39_repeat_masker.gtf"
    build_index "${REF_MUS}/Mus_musculus.GRCm39.dna.primary_assembly.fa" \
                "${REF_MUS}/Mus_musculus.GRCm39.115.gtf" "${REF_MUS}/star_index_mus"
fi

# ===========================================================================
# 3. STAGE 2 -- SRA RUN TABLE + GSM->SRR MAP  (same logic as the Snakefile)
# ===========================================================================
echo "----- STAGE 2: metadata -----"
RUN_TABLE="${BASE_DIR}/SraRunTable.csv"
MAP_FILE="${BASE_DIR}/gsm_srr_map.tsv"     # lines: GSM<TAB>SRR1,SRR2,...

python - "$SRP_ID" "$RUN_TABLE" "$MAP_FILE" <<'PY'
import sys, os, json, csv, urllib.request, urllib.parse
srp, run_table, map_file = sys.argv[1], sys.argv[2], sys.argv[3]
if not os.path.exists(run_table):
    print(f"Fetching run table for {srp} from NCBI ...")
    sp = {"db": "sra", "term": srp, "retmode": "json", "retmax": 500}
    su = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(sp)
    ids = json.loads(urllib.request.urlopen(su).read().decode()).get("esearchresult", {}).get("idlist", [])
    if ids:
        fp = {"db": "sra", "id": ",".join(ids), "rettype": "runinfo", "retmode": "text"}
        fu = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(fp)
        urllib.request.urlretrieve(fu, run_table)
gsm = {}
with open(run_table) as f:
    delim = '\t' if '\t' in f.readline() else ','
with open(run_table) as f:
    r = csv.DictReader(f, delimiter=delim); h = r.fieldnames
    rc = next((c for c in h if c.lower() in ('run', 'run_s')), None)
    gc = next((c for c in h if c.lower() == 'samplename'), None)
    if not (rc and gc):
        sys.exit("ERROR: could not find Run / SampleName columns in run table")
    for row in r:
        s, g = row[rc].strip(), row[gc].strip()
        if s and g:
            gsm.setdefault(g, []).append(s)
if not gsm:
    sys.exit("ERROR: no GSM/SRR pairs parsed (NCBI fetch failed?)")
with open(map_file, 'w') as o:
    for g in sorted(gsm):
        o.write(g + '\t' + ','.join(gsm[g]) + '\n')
print(f"Mapped {len(gsm)} GSMs / {sum(len(v) for v in gsm.values())} SRRs")
PY

cut -f1 "$MAP_FILE" > "${BASE_DIR}/GSM_Acc_List.txt"

# ===========================================================================
# 4. STAGE 3 -- PER-SAMPLE: download -> STAR Solo -> velocyto
# ===========================================================================
# To process several samples concurrently on one fat node instead of serially,
# replace the `while` loop with GNU parallel, e.g.:
#     export -f process_gsm ; ... ; parallel -j 2 ...
# but keep -j small: each velocyto needs ~160 GB RAM with these defaults.
echo "----- STAGE 3: per-sample processing -----"
STAR_DIR="${BASE_DIR}/${SRP_ID}_STAR"
VELO_DIR="${BASE_DIR}/${SRP_ID}_velocyto"
LOOMS=()

while IFS=$'\t' read -r gsm srrs; do
    loom="${VELO_DIR}/${gsm}_velocyto/${gsm}_Aligned.sortedByCoord.out.loom"
    LOOMS+=("$loom")

    if [ -f "$loom" ]; then
        echo "  [skip] ${gsm}: loom already exists"
        continue
    fi

    bam="${STAR_DIR}/${gsm}_Aligned.sortedByCoord.out.bam"
    barcodes="${STAR_DIR}/${gsm}_Solo.out/Gene/filtered/barcodes.tsv"

    # ---- download + concatenate + align (skip if BAM already present) -------
    if [ ! -f "$bam" ]; then
        tmp="${BASE_DIR}/fastq/${gsm}_temp"
        r1="${BASE_DIR}/combined_fastq/${gsm}_1.fastq.gz"
        r2="${BASE_DIR}/combined_fastq/${gsm}_2.fastq.gz"
        mkdir -p "$tmp" "$(dirname "$r1")"
        : > "$r1"; : > "$r2"

        echo "  [${gsm}] downloading runs: ${srrs}"
        IFS=',' read -ra arr <<< "$srrs"
        for srr in "${arr[@]}"; do
            prefetch --max-size 100G -O "$tmp" "$srr"
            fasterq-dump -e "$DL_THREADS" --outdir "$tmp" "$tmp/$srr"
            pigz -f -p "$DL_THREADS" "$tmp/${srr}_1.fastq"
            pigz -f -p "$DL_THREADS" "$tmp/${srr}_2.fastq"
            cat "$tmp/${srr}_1.fastq.gz" >> "$r1"
            cat "$tmp/${srr}_2.fastq.gz" >> "$r2"
            rm -rf "$tmp/$srr" "$tmp/${srr}_1.fastq.gz" "$tmp/${srr}_2.fastq.gz"
        done
        rm -rf "$tmp"

        echo "  [${gsm}] STAR Solo alignment"
        mkdir -p "$STAR_DIR"
        STAR --runThreadN "$STAR_THREADS" \
             --genomeDir "$GENOME_DIR" \
             --readFilesCommand zcat \
             --outFileNamePrefix "${STAR_DIR}/${gsm}_" \
             --readFilesIn "$r2" "$r1" \
             --soloType CB_UMI_Simple \
             --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
             --soloBarcodeReadLength 0 \
             --soloCBwhitelist "$WHITELIST_V2" \
             --soloCellFilter EmptyDrops_CR \
             --soloStrand Forward \
             --outSAMattributes NH HI AS nM CB UB \
             --outSAMtype BAM SortedByCoordinate

        # combined fastqs were Snakemake temp() -> drop them to save space
        rm -f "$r1" "$r2"
    else
        echo "  [${gsm}] BAM exists -> skipping download + STAR"
    fi

    # ---- velocyto ----------------------------------------------------------
    echo "  [${gsm}] velocyto run"
    out_dir="${VELO_DIR}/${gsm}_velocyto"
    mkdir -p "$out_dir"
    velocyto run -b "$barcodes" \
                 -o "$out_dir" \
                 -m "$MASK_FILE" \
                 -e "${gsm}_Aligned.sortedByCoord.out" \
                 -@ "$VELO_THREADS" \
                 --samtools-memory 40000 \
                 "$bam" \
                 "$GTF_FILE"
done < "$MAP_FILE"

# ===========================================================================
# 5. STAGE 4 -- AGGREGATION + QC
# ===========================================================================
echo "----- STAGE 4: aggregation + QC -----"
# Resolve the CXXABI issue the same way the Snakefile does (env is active here)
export LD_PRELOAD="${CONDA_PREFIX}/lib/libstdc++.so.6"

python data_pp_smk.py \
    --inputs "${LOOMS[@]}" \
    --output "${BASE_DIR}/${SRP_ID}_processed.h5ad" \
    --srp "$SRP_ID" \
    --qc_dir "$PP_QC_DIR"

echo "=========================================================="
echo "  DONE -> ${BASE_DIR}/${SRP_ID}_processed.h5ad"
echo "  $(date)"
echo "=========================================================="
