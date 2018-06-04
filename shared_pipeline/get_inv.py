import numpy as np

data0=open('hg38.200bp.noblacklist_id.bed','r')
data01={}
data02=[]

i=0
for records in data0:
	info = records.split()
	if not (info[0] in data01):
		data01[info[0]] = i
		data02.append(info[0])
	i=i+1
	if i% 100000 ==0:
		print(i)



data0.close()


info_table = []
for chrom in data02:
	info = data01[chrom]
	print([chrom, info])
	info_table.append([chrom, info])


info_table = np.array(info_table)

end = info_table[1:,1]
end = np.append(end, i-1)
end = end.reshape(end.shape[0], 1)

info_table = np.concatenate((info_table, end), 1)

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()


write2d_array(info_table,'hg38.200bp.noblacklist_id.inv' )

