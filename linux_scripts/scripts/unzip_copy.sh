# The script is used to unzip the .gz format files and generate a new copy from the origin fastq.gz files

find /data/rsgeno1/liuzhen -name "*.fastq.gz" | while read line
do
	path=${line%/*}
	filename=${line##*/}
	new_path="/data/rsgeno1/liuzhen/new_Caucasian${path:35}"
	new_filename=${filename%.gz}
	new_file="$new_path/$new_filename"
	mkdir -p $new_path
	gzip -cd $line > $new_file
	printf "Finish gunzip $new_file"
done
