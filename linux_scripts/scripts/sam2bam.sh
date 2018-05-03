# The script is used to transfer sam files to bam files

for sam_file in $(find /data/rsgeno1/liuzhen/new_Caucasian -name "*.sam")
do
	samtools view -bS -o "${sam_file%.*}.bam" $sam_file
done
