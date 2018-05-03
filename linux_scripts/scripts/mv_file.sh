# The script is used to move a kind of files from one directory to another

find /data/rsgeno1/liuzhen/new_Caucasian -name "*_R*" | while read dir_path
do
	mod_type="${dir_path#*adapter_quality_trimmed/}"
	new_dir_path="${dir_path%/adapter_quality_trimmed*}/adapter_quality_trimmed/${mod_type%/*}_${mod_type#*/}"
	mv $dir_path $new_dir_path
done
