# The script is used to convert the _treat_pileup.bdg files to bigwig files. The _treat_pileup.bdg files is generate by peak calling of MACS2. The bdg files will be sorted first and the reference genome is hg19

for ChIP_dir in $(find /data/rsgeno1/liuzhen/new_Caucasian -name "H3K*")
do
        for SRR_dir in $(find "${ChIP_dir}/call_peak" -name "*_R*")
        do
		for bdg_file in $(find "${SRR_dir}" -name "*.bdg")
		do
			bdg_file_name=${bdg_file%.bdg}
			printf "Sorting ${bdg_file}\n"
			LC_COLLATE=C sort -k1,1 -k2,2n $bdg_file > "${bdg_file_name}.sort.bdg"
			printf "Converting ${bdg_file_name}.sort.bdg to bigwig file\n"
			bedGraphToBigWig "${bdg_file_name}.sort.bdg" /data/rsgeno1/liuzhen/genome/hg19_chrSize.txt "${bdg_file_name}.bw"
			printf "Finish converting\n"
		done
	done
done

