# The script is used to process bam files for calling peaks

control="/data/rsgeno1/liuzhen/new_Caucasian/GM12891/Input/mapping_sam/GM12891_Input_R1/SRR998238.bam"
for ChIP_dir in $(find /data/rsgeno1/liuzhen/new_Caucasian/GM12891 -name "H3K*")
do
	for SRR_dir in $(find "${ChIP_dir}/mapping_sam" -name "*_R*")
	do
		new_dir="${ChIP_dir}/call_peak/${SRR_dir##*mapping_sam/}"
		mkdir -p $new_dir
		bam_file=$(find $SRR_dir -name "*.bam")
		tmp_name=${bam_file##*/}
		SRR_num=${tmp_name%.bam}
		printf "MACS2 callpeak running ${bam_file}......\n"
		printf "Control is ${control}\n"
		macs2 callpeak -t $bam_file -c $control -g hs -n $SRR_num -B -q 0.01 --outdir $new_dir > ${new_dir}/${SRR_num}_MACS2_log 2>&1
		printf "MACS2 callpeak finish."
	done
done
