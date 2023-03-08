#!/bin/bash

# Declare vars. 
rna_bam=""
patient_id=""
rna_sample_id=""
mut_path=""
out_dir=""
work_dir=""

while getopts "b:p:m:o:w:s" opt; do
    case "${opt}" in
        b)
            rna_bam="${OPTARG}";
            ;;
        p)
            patient_id="${OPTARG}";
            ;;   
        m)
            mut_path="${OPTARG}";
            ;;                         
        o)
            out_dir="${OPTARG}";
            ;;
        w)
            work_dir="${OPTARG}";
            ;;
        s)
            rna_sample_id="${OPTARG}";
            ;;                        
        \?)
            echo -e "Error: Invalid option: -${OPTARG}\n" >&2
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))

# Check if required arguments are set
if [ -z "$rna_bam" ] || \
   [ -z "$patient_id" ] || \
   [ -z "$rna_sample_id" ] || \
   [ -z "$mut_path" ] || \
   [ -z "$out_dir" ] || \
   [ -z "$work_dir" ]
then
    echo -e "Error: Missing required argument(s)\n"
    echo -e "Usage: ./get_RNA_VAF.sh -b <rna_bam> -p <patient_id> -s <rna_sample_id> -m <mut_path> -o <out_dir> -w <work_dir>\n"
    exit 1
fi

# paths
baseDir="${PWD}"
dbsnp_138_hg38="${baseDir}/assets/Homo_sapiens_assembly38.dbsnp138.vcf"
singularityDir="${baseDir}/singularity_images/"
hg38_reference="${baseDir}/assets/GRCh38.d1.vd1.fa"

singularity_command="singularity exec --bind ${PWD}:${PWD} ${singularityDir}/rvdriver_latest.sif"

# Check if the working directory exists
if [ ! -d "$work_dir" ]; then
  # If the directory doesn't exist, create it
  mkdir -p "$work_dir"
fi

# Check if the out directory exists
if [ ! -d "$out_dir" ]; then
  # If the directory doesn't exist, create it
  mkdir -p "$out_dir"
fi

# Change into the directory
cd "$work_dir"

# make tmp_dir
mkdir ./tmp_dir

echo working on ${patient_id}

# check if paired end or single end
paired_test=$( { samtools view -H "${rna_bam}" ; samtools view "${rna_bam}" | head -n1000; } | samtools view -c -f 1)

echo "first 1000 lines of rna_bam contain ${paired_test} paired reads"

if [ "$paired_test" -eq "0" ]; then
    echo "${rna_bam} is single end"
    # keep only mapped reads
    echo "filtering bam for proper mapped reads"
    ${singularity_command} samtools view -b -F 0x4 "${rna_bam}" > "${rna_sample_id}"_mapped_pairs.bam
    echo "DONE!"
else
    echo "${rna_bam} is paired end"
    # Keep only reads in proper mapped pair
    echo "filtering bam for proper mapped pairs"
    ${singularity_command} samtools view -b -f 2 "${rna_bam}" > "${rna_sample_id}"_mapped_pairs.bam
    echo "DONE!"
fi

#Run the first step: Mark duplicates
${singularity_command} gatk --java-options '-Djava.io.tmpdir=./tmp_dir -Xmx8g -Xms4g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1' MarkDuplicates \
    -I "${rna_sample_id}"_mapped_pairs.bam \
    -O "${rna_sample_id}"_md.bam \
    -M "${rna_sample_id}"_md_metrics.txt \
    --TMP_DIR ./tmp_dir

#Second step: Split N Cigar reads. Split reads with Ns over introns
${singularity_command} gatk --java-options '-Djava.io.tmpdir=./tmp_dir -Xmx8g -Xms4g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1' SplitNCigarReads \
    -R "${hg38_reference}" \
    -I "${rna_sample_id}"_md.bam \
    -O "${rna_sample_id}"_md_split.bam \
    --tmp-dir ./tmp_dir

#Third step: Base recalibration. Detect and correct for patterns of systematic errors in the base quality scores.
#Generate recalibration stats
${singularity_command} gatk --java-options '-Djava.io.tmpdir=./tmp_dir -Xmx8g -Xms4g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1' BaseRecalibrator \
    -I "${rna_sample_id}"_md_split.bam \
    -R "${hg38_reference}" \
    --known-sites "${dbsnp_138_hg38}" \
    -O "${rna_sample_id}".table \
    --tmp-dir ./tmp_dir

#Apply recalibration stats
${singularity_command} gatk --java-options '-Djava.io.tmpdir=./tmp_dir -Xmx8g -Xms4g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1' ApplyBQSR \
    -R "${hg38_reference}" \
    -I "${rna_sample_id}"_md_split.bam \
    --bqsr-recal-file "${rna_sample_id}".table \
    -O pre_processed_"${rna_sample_id}".bam \
    --tmp-dir ./tmp_dir

# Remove the old index file and tmp directory
rm -f pre_processed_"${rna_sample_id}"*bai
rm -rf ./tmp_dir

# Index the new BAM file
echo "Indexing preprocessed bam"
${singularity_command} samtools index pre_processed_"${rna_sample_id}".bam
echo "DONE!"

rna_bam="pre_processed_${rna_sample_id}.bam"

# Remove intermeidate files
ls "${rna_sample_id}"_md* | xargs rm
rm "${rna_sample_id}"_mapped_pairs.bam

# get mutation expression 

${singularity_command} Rscript "${baseDir}"/bin/mut_expression.R \
    --patient_id "${patient_id}" --rna_bam "${rna_bam}" \
    --sample_id "${rna_sample_id}" --out_dir "${out_dir}" --mut_path "${mutation_table}"
