# The script is used to convert .bam files to .bigwig files using deeptools. The unsorted .bam files  will be sorted by samtools and .bai index files will generate.

for RNAseq_dir in $(find /data/rsgeno1/liuzhen/new_Caucasian -name "RNA*")
do
        for SRR_dir in $(find "${RNAseq_dir}/STAR_mapping" -name "*_R*")
        do
                for bam_file in $(find "${SRR_dir}" -name "*.bam")
                do
                        bam_file_name=${bam_file%.bam}
                        printf "Samtools sorting ${bam_file}\n"
                        samtools sort -m 16G $bam_file -o "${bam_file_name}.sort.bam"
                        samtools index "${bam_file_name}.sort.bam"
                        printf "Converting ${bam_file_name}.sort.bam to bigwig file\n"
                        bamCoverage -b "${bam_file_name}.sort.bam" -o "${bam_file_name}.bw"
                        printf "Finish converting\n"
                        printf "\n"
                done
        done
done