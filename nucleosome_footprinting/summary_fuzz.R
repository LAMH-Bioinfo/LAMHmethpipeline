library(data.table)
df_sample <- read.table('../sample_list_2.txt')
list_sample <- as.vector(unlist(df_sample$V1))
batch_id <- 'Colon_case_control_780samples'
for(i in 1:length(list_sample)) {
df <- fread(paste0(list_sample[i], '.outMean.fuzziness.tsv'), header=TRUE)
df <- data.frame(df)
df <- df[,c(1,2)]
colnames(df) <- c('cpg_id', list_sample[i])
if(i==1) {
df_out <- df
} else {
df_out <- merge(df_out, df, by='cpg_id', all=TRUE)
}
}
write.table(df_out, paste0(batch_id, '_fuzziness.tsv'), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')



