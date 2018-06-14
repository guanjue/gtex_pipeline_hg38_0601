input_folder=/storage/home/gzx103/scratch/gtex_encode/bams/
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
#ref0=$(head -1 all_ref_list.txt | awk -F '.' -v OFS='\t' '{print $1"."$2"."$3"."$4".nbp_2r_bgadj.txt"}')
ref0=$(head -1 all_ref_list.txt | awk -F '_' -v OFS='_' '{print $1,$2,$3,$4".fisher_p.txt"}')

###### (5) PKnorm normalization (between reference sample)
while read LINE
do
	sig1=$ref0
	#sig2=$(echo "$LINE" | awk -F '.' -v OFS='\t' '{print $1"."$2"."$3"."$4".nbp_2r_bgadj.txt"}')
	sig2=$(echo "$LINE" | awk -F '_' -v OFS='_' '{print $1,$2,$3,$4".fisher_p.txt"}')
	sig2_celltype=$(echo "$LINE" | awk -F '_' -v OFS='_' '{print $1,$2,$3,$4}')
	upperlim=1000
	lowerlim=0
	echo $sig1 
	echo $sig2
	echo $sig2_celltype
	### set upper limit
	cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
	cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt'
	### peak norm
	time python $script_folder'peaknorm_rotate_log_ref_mean.py' -w 1 -p 1 -n 500000 -a 1 -b $sig1'.upperlim.txt' -c 1 -d $sig2'.upperlim.txt' -u $upperlim -l $lowerlim
	### rm tmp files
	rm $sig1'.upperlim.txt' $sig2'.upperlim.txt'
	mv $sig2_celltype'.pknorm.txt' $sig2_celltype'.pknorm.ref.txt'
	mv $sig2_celltype'.info.txt' $sig2_celltype'.info.ref.txt'
	mv $sig2_celltype'.pknorm.scatterplot.png' $sig2_celltype'.pknorm.scatterplot.ref.png'
	mv $sig2_celltype'.scatterplot.png' $sig2_celltype'.scatterplot.ref.png'
done < all_ref_list.txt
ref0_celltype=$(echo "$ref0" | awk -F '.' -v OFS='\t' '{print $1}' | awk -F '_' -v OFS='_' '{print $1,$2,$3,$4}')
cp $ref0 $ref0_celltype'.pknorm.ref.txt'


###### (6) PKnorm normalization (between reference & target sample)
for mk in $(cat mark_list.txt)
do
	echo $mk
	ref_mk=$(head -1 $mk'.ref_frip.txt' | awk -F '_' -v OFS='_' '{print $1,$2,$3,$4".fisher_p.ref.txt"}')
	cat $mk'.file_list.txt' | awk -F '_' -v OFS='_' '{print $1,$2,$3,$4}' > $mk'.file_list.ct.txt'
	echo $ref_mk
	ls $ref_mk
	while read LINE
	do
		sig1=$ref_mk
		sig2=$LINE'.fisher_p.txt'
		sig2_celltype=$(echo "$LINE" | awk -F '.' -v OFS='\t' '{print $1}')
		upperlim=1000
		lowerlim=0
		echo $sig1 
		echo $sig2
		echo $sig2_celltype
		### set upper limit
		cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
		cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt'
		### peak norm
		time python $script_folder'peaknorm_rotate_log_z_mean.py' -w 1 -p 1 -n 500000 -a 1 -b $sig1'.upperlim.txt' -c 1 -d $sig2'.upperlim.txt' -u $upperlim -l $lowerlim
		### rm tmp files
		rm $sig1'.upperlim.txt' $sig2'.upperlim.txt'
		cat $sig2_celltype'.pknorm.txt' | awk -F '\t' -v OFS='\t' '{if ($1>=16) print 16; else if ($1<=2) print 2; else print $1}' > $sig2_celltype'.pknorm.2_16lim.txt'
	done < $mk'.file_list.ct.txt'
done

#mkdir entex_data_output_2_16lim/
#mv *.pknorm.2_16lim.txt entex_data_output_2_16lim/

mkdir entex_data_output_2_16lim_fisherp/
mv *.pknorm.2_16lim.txt entex_data_output_2_16lim_fisherp/

mkdir entex_data_output_0_16lim_fisherp/
mv *.pknorm.0_16lim.txt entex_data_output_0_16lim_fisherp/

mkdir pknorm_ref_fisherp
mv *.pknorm.ref.txt pknorm_ref_fisherp/
mv *.info.ref.txt pknorm_ref_fisherp/
mv *ref.png pknorm_ref_fisherp/

mkdir pknorm_fisherp
mv *png pknorm_fisherp/
mv *.pknorm.txt pknorm_fisherp/
mv *.info.txt pknorm_fisherp/

mkdir fisher_p
mv *.fisher_p.txt fisher_p/

mkdir nbp_2r
mv *.signal.tab.nbp_2r_bgadj.txt nbp_2r/
mv *.signal.tab.mvsp.txt nbp_2r/

mkdir reads_count
mv *.bam.signal.tab reads_count/
mv *.bam.tab reads_count/

###### write ideas input file
ls entex_data_output_2_16lim_fisherp/ > entex_data_output_2_16lim_fisherp_filelist.txt
### remove previous ideas.input file
if [ -f ideas.2_16lim_fisherp.input ]
then
	echo 'remove previous ideas.input'
	rm ideas.2_16lim_fisherp.input
fi
### get new ideas.input file
for filename in $(cat entex_data_output_2_16lim_fisherp_filelist.txt)
do
	echo $filename
	mk=$(echo "$filename" | awk -F '_' '{print $2}' | awk -F '-' '{print $1}')
	exp=$(echo "$filename" | awk -F '_' '{print $1}' | awk -F '-' '{print $1}')
	if [ $exp = "DNase" ]; then mk=$(echo "$filename" | awk -F '_' '{print $1}' | awk -F '-' '{print $1}');fi
	ct=$(echo "$filename" | awk -F '_' '{print $4}' | awk -F '.' '{print $1}' )
	rep=$(echo "$filename" | awk -F '.' '{print $1}' | awk -F '_' '{print $3}')
	echo $ct $mk $input_folder'entex_data_output_2_16lim_fisherp/'$filename >> ideas.2_16lim_fisherp.input
done


###### write ideas input file
ls entex_data_output_0_16lim_fisherp/ > entex_data_output_0_16lim_fisherp_filelist.txt
### remove previous ideas.input file
if [ -f ideas.0_16lim_fisherp.input ]
then
	echo 'remove previous ideas.input'
	rm ideas.0_16lim_fisherp.input
fi
### get new ideas.input file
for filename in $(cat entex_data_output_0_16lim_fisherp_filelist.txt)
do
	echo $filename
	mk=$(echo "$filename" | awk -F '_' '{print $2}' | awk -F '-' '{print $1}')
	exp=$(echo "$filename" | awk -F '_' '{print $1}' | awk -F '-' '{print $1}')
	if [ $exp = "DNase" ]; then mk=$(echo "$filename" | awk -F '_' '{print $1}' | awk -F '-' '{print $1}');fi
	ct=$(echo "$filename" | awk -F '_' '{print $4}' | awk -F '.' '{print $1}' )
	rep=$(echo "$filename" | awk -F '.' '{print $1}' | awk -F '_' '{print $3}')
	echo $ct $mk $input_folder'entex_data_output_0_16lim_fisherp/'$filename >> ideas.0_16lim_fisherp.input
done


###### write ideas input file
ls pknorm_fisherp/*pknorm.txt > entex_data_output_0_1000_log2_fisherp_filelist.txt
### remove previous ideas.input file
if [ -f ideas.0_1000_log2_fisherp.input ]
then
	echo 'remove previous ideas.input'
	rm ideas.0_1000_log2_fisherp.input
fi
### get new ideas.input file
for filename in $(cat entex_data_output_0_1000_log2_fisherp_filelist.txt)
do
	echo $filename
	mk=$(echo "$filename" | awk -F '_' '{print $2}' | awk -F '-' '{print $1}')
	exp=$(echo "$filename" | awk -F '_' '{print $1}' | awk -F '-' '{print $1}')
	if [ $exp = "DNase" ]; then mk=$(echo "$filename" | awk -F '_' '{print $1}' | awk -F '-' '{print $1}');fi
	ct=$(echo "$filename" | awk -F '_' '{print $4}' | awk -F '.' '{print $1}' )
	rep=$(echo "$filename" | awk -F '.' '{print $1}' | awk -F '_' '{print $3}')
	echo $ct $mk $input_folder'pknorm_fisherp/'$filename >> ideas.0_1000_log2_fisherp.input
done



paste cell_list.txt cell_list.txt cell_list.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,"255,255,255"}' > cellinfo.txt

mkdir entex_data_output_2_16lim_ideas

###### run IDEAS
time Rscript bin/runme.R ideas.input entex.parafile entex_data_output_2_16lim_ideas

###bed
#row_number: 15441243
for file in $(cat file_list.txt)
do
	echo $file
	head -15441243 $file > $file'.tmp'
	mv $file'.tmp' $file
done

