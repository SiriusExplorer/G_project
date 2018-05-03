# The script is used to map RNA-seq on hg19 genome using STAR

for RNA_seq_dir in $(find /data/rsgeno1/liuzhen/new_Caucasian -name "*RNA_seq")
do
	for SRR_dir in $(find $RNA_seq_dir -name "*_R*")
	do
		new_dir="${SRR_dir%/adapter*}/STAR_mapping/${SRR_dir##*trimmed/}"
                mkdir -p $new_dir
		printf "STAR mapping ${SRR_dir##*trimmed/}......\n"
		STAR --runThreadN 20 --genomeDir /data/rsgeno1/liuzhen/index/hg19_star --readFilesIn $SRR_dir/*.fq --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$new_dir/"
		printf "STAR mapping ${SRR_dir##*trimmed/} finnish\n"
	done
done
