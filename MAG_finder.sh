#!/bin/bash

# Function to display help
display_help() {
    echo "Usage: $0 -i <input_directory> -o <output_directory> -t <number_of_threads>"
    echo
    echo "This script performs contig binning, check results and find a taxon representative automatically"
    echo "MAG_finder use several different binning methods and softwares"
    echo "Bin results are saved to the output directory."
    echo
    echo "Arguments:"
    echo "  -i   Path to the directory containing read files (_1.fq and _2.fq)."
    echo "  -o   Path to the directory where results will be saved."
    echo "  -t   Number of threads to be used by the tools."
    echo
}

# Function for logging
log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1" | tee -a "$log_file"
}

# Check if no argument was provided
if [ $# -eq 0 ]; then
    display_help
    exit 1
fi

# Parse command-line arguments
while getopts i:o:t:h opt; do
  case $opt in
    i) reads_dir=$OPTARG ;;
    o) out_dir=$OPTARG ;;
    t) threads=$OPTARG ;;
    h) display_help
       exit 0 ;;
    *) display_help
       exit 1 ;;
  esac
done

# Obtenha o caminho do diretório onde o script está sendo executado
script_dir=$(dirname "$(dirname "$(readlink -f "$0")")")

# Verify if required arguments were provided
if [ -z "$reads_dir" ] || [ -z "$out_dir" ] || [ -z "$threads" ]; then
    log "Missing arguments. See usage below:"
    display_help
    exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
    log "Created output directory: $out_dir"
fi

# Create log file
log_file="${out_dir}/log.txt"
touch "$log_file"

# FastQC quality check
log "Starting FastQC quality check..."
mkdir -p "${out_dir}/01_FASTQC"
fastqc -t $threads -o "${out_dir}/01_FASTQC" "${reads_dir}"/*

# Trimmomatic for quality trimming
log "Starting Trimmomatic for quality trimming..."
mkdir -p "${out_dir}/02_TRIMMED_READS"
for reads_1 in "${reads_dir}"/*_1.*; do
    sample=$(basename "$reads_1")
    sample=${sample%_1.*}
    read_file_1=$reads_1
    read_file_2=${reads_1/_1./_2.}
    trimmed_file_1="${out_dir}/02_TRIMMED_READS/${sample}_trimmed_1.fastq"
    trimmed_file_2="${out_dir}/02_TRIMMED_READS/${sample}_trimmed_2.fastq"
    trimmomatic PE -threads $threads $read_file_1 $read_file_2 $trimmed_file_1 "${out_dir}/02_TRIMMED_READS/${sample}_trimmed_1_unpaired.fastq" $trimmed_file_2 "${out_dir}/02_TRIMMED_READS/${sample}_trimmed_2_unpaired.fastq" SLIDINGWINDOW:4:20 MINLEN:36
done

# MEGAHIT for assembly
log "Starting MEGAHIT for assembly..."
mkdir -p "${out_dir}/03_ASSEMBLY"
for trimmed_file_1 in "${out_dir}/02_TRIMMED_READS/"*_trimmed_1.fastq; do
    sample=$(basename "$trimmed_file_1")
    sample=${sample%_trimmed_1.fastq}
    trimmed_file_2="${trimmed_file_1/_trimmed_1.fq/_trimmed_2.fastq}"
    megahit -1 $trimmed_file_1 -2 $trimmed_file_2 -t $threads -o "${out_dir}/03_ASSEMBLY/${sample}"
done

log "MAG finding process completed. Check ${out_dir} for results."

# calculating coverage
mkdir -p "${out_dir}/04_COVERAGE"
bash $script_dir/MAG_finder/gen_cov_file.sh -a ${out_dir}/03_ASSEMBLY/${sample}/*.fa -o ${out_dir}/04_COVERAGE/ "${out_dir}/02_TRIMMED_READS/"*_trimmed_{1,2}.fastq -t $threads -m 80 

#Metabat2
for i in ${out_dir}/04_COVERAGE/work_files/*.bam; do
	$name=$(basename "$i" | cut -d '_' -f 1) 
	jgi_summarize_bam_contig_depths --outputDepth ${out_dir}/04_COVERAGE/work_files/$name.txt $i; 
	metabat2 -i ${out_dir}/03_ASSEMBLY/${name}.fa -o bin -a ${out_dir}/04_COVERAGE/work_files/${name}.txt
done



