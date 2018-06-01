###### (1) Reads-count --> (2) NB-p-value 
for ct in $(cat cell_list.txt)
do
	echo $ct
	for mk in $(cat mark_list.txt)
	do
		echo $mk
		ls *$mk*$ct*.signal.tab > file_list_tmp.txt
		###### 2_nbp
		for file in $(cat file_list_tmp.txt)
		do
			echo $file
			time Rscript $script_folder'negative_binomial_p_2r_bgadj_frip.R' $file $file
		done
	done
done
rm file_list_tmp.txt


###### extract all reference list
for mk in $(cat mark_list.txt)
do
	echo $mk
	ls *$mk*.signal.tab.mvsp.txt > $mk'.file_list.txt'
	time Rscript $script_folder'get_mk_ref.R' $mk'.file_list.txt' $mk'.ref_frip.txt'
done

###### get all reference list
cat *.ref_frip.txt | sort -k2,2rn > all_ref_list.txt
ref0=$(head -1 all_ref_list.txt | cut -f2)


######






file=DNase-seq__bam_vagina_ENCSR437AYW_1_ENCFF244BRN.bam
time ~/group/software/ucsc/bigWigAverageOverBed $file'.bw' /storage/home/gzx103/scratch/gtex_encode/hg38.200bp.noblacklist_id.bed $file'.tab'
cut -f5 $file'.tab' > $file'.signal.tab'


