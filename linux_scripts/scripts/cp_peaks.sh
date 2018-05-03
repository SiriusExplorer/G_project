# The script is used to copy a kind of files in a specified dirctory

for xls_file in $(find /data/rsgeno1/liuzhen/new_Caucasian -name "*xls")
do
	xls_file_prefix=${xls_file%/*}
	xls_file_name=${xls_file_prefix##*/}
	cp $xls_file "/home/liuzhen/G_project/data/MACS2_peaks/${xls_file_name}.xls"
done
