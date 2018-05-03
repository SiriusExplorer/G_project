# The script is used to run bowtie and map .fq to the hg19 index

find /data/rsgeno1/liuzhen/new_Caucasian -name "adapter*" -type d | while read trim_dir
do
	find $trim_dir -name "*_R*" -type d | while read SRR_dir
	do
		new_dir="${SRR_dir%/adapter*}/mapping_sam/${SRR_dir##*trimmed/}"
		mkdir -p $new_dir
		flag=0
		find $SRR_dir -name "*.fq" > ~/tmp/find_tmp
		while read fq_file
		do
			tmp=${fq_file##*/}
			SRR_num=${tmp%%_*}
			if [ $flag == 0 ]
			then
				fq_file1=$fq_file
				old_SRR_num=$SRR_num
			else
				fq_file2=$fq_file
			fi
			let "flag++"
			if [ $old_SRR_num != $SRR_num ]
			then
				printf "Error! The SRR number is different in one SRR directory\n"
			fi
		done < ~/tmp/find_tmp
		printf "Bowtie mapping ${SRR_dir} ${SRR_num}......\nThe fq_file is:\n${fq_file1}\n${fq_file2}\n"
		bowtie -n 2 --nomaqround -I 50 -X 550 --fr --maxbts 200 -k 1 -m 1 --chunkmbs 512 -p 15 -t /data/rsgeno1/liuzhen/index/hg19/hg19 -1 $fq_file1 -2 $fq_file2 -S "${new_dir}/${SRR_num}.sam" > "${new_dir}/${SRR_num}_mapping_log" 2>&1 
		printf "Bowtie mapping finish ${SRR_num}\n."
	done
done
