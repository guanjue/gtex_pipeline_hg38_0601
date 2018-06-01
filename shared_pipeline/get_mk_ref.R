### get parameters
args = commandArgs(trailingOnly=TRUE)

parameters_file_list = args[1]
output_name = args[2]

parameters_files = read.table(parameters_file_list, header=F)

FRiP_list = c()
i = 0
for (file in parameters_files){
	i = i+1
	print(i)
	parameters = read.table(file, header=F)
	### FRiP score is in the 5th column
	FRiP_list[i] = parameters[5]
}


ref_file = parameters_files[which.max(FRiP_list)]
frip_ref = FRiP_list[which.max(FRiP_list)]

write.table(c(output_name, frip_ref), sep='\t', quote=F, col.names=F, row.names=F)
