# The script is used to run the fastqc analyze to check the quality of the fastq file

find /data/rsgeno1/liuzhen/new_Caucasian -name "*.fastq" | while read line
do
        path=${line%/*}
	fastqc -o $path $line
	printf "Finish FastQC $line"
done
