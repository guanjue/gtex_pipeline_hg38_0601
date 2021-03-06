### get parameters
args = commandArgs(trailingOnly=TRUE)

parameters_file_list = args[1]
output_name = args[2]

parameters_files = read.table(parameters_file_list, header=F)

FRiP_list = c()
for (i in c(1:dim(parameters_files)[1])){
	print(i)
	print(parameters_files[i,1])
	file=toString(parameters_files[i,1])
	parameters = read.table(file, header=F)
	### FRiP score is in the 5th column
	FRiP_list[i] = parameters[5]
}


ref_file = toString(parameters_files[which.max(FRiP_list),1])
frip_ref = FRiP_list[which.max(FRiP_list)]

write.table(c(ref_file, frip_ref), output_name, sep='\t', quote=F, col.names=F, row.names=F)
