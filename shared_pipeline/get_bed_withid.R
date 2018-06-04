file = '/storage/home/gzx103/scratch/gtex_encode/hg38.200bp.noblacklist.bed'
output = '/storage/home/gzx103/scratch/gtex_encode/hg38.200bp.noblacklist_id.bed'

data = read.table(file, header=F)

num = dim(data)[1]

id = c(1:num)

data_id = cbind(data, id)

write.table(data_id, output, sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
