# The script is used to zip the files

find /data/rsgeno1/liuzhen/new_Caucasian -name "*.fastq" | while read line
do
	gzip $line
	printf "Finish gzip $line"
done
