# The script is used to run trim_galore to cutadapter of the fastq files and execute quality control

n=1
sample_mod="Initial_value"
find /data/rsgeno1/liuzhen/new_Caucasian -name "SRR*" -type d | while read SRR_path
do
	old_sample_mod=$sample_mod
	sample_mod_temp=${SRR_path%/raw_data*}
	sample_mod=${sample_mod_temp##*Caucasian/}
	if [ $sample_mod != $old_sample_mod ]
	then
		n=1
	fi
	new_path="${SRR_path%/raw_data*}/adapter_quality_trimmed/${sample_mod%/*}_${sample_mod#*/}_R${n}"
	mkdir -p $new_path
	echo "Trim_galore Processing ${SRR_path##*/}..."
	trim_galore --fastqc -o $new_path --trim1 --paired $(find $SRR_path -name "*.fastq") > "${new_path}/${SRR_path##*/}_trimmng_log" 2>&1
	echo "Trim_galore Finish ${SRR_path##*/}"
	let "n++"
done
