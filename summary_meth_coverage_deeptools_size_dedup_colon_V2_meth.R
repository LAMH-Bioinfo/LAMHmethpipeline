library(data.table)
sample_list <- read.table('../sample_id.txt', header = FALSE)
sample_list <- sample_list$V1
sample_name <- 'Standards'
work_dir <- '/shared/MingmaData/Project_s1867g02012_21Samples_20220213_1644737104/all_meth_deeptools_80_200/'
file_name_meth <- paste(sample_name, '_methratio_deeptools_size_80_200_sorted.csv', sep = '')
file_name_coverage <- paste(sample_name, '_coverage_deeptools_size_80_200_sorted.csv', sep = '')
file_appendix <- '_meth_combined_length_80_200_selected.bed'

for(sample in sample_list) {
  df <- read.table(paste(work_dir, sample, file_appendix, sep = ''), header = FALSE)
  df <- df[,c(1,2,3,4,5,6,7)]
  colnames(df) <- c('chr_probe', 'start_probe', 'end_probe', 'chr', 'start', 'end', sample)
  df$probe_id <- paste(df$chr_probe, df$start_probe, df$end_probe, sep='_')
  df <- df[,c('probe_id','chr_probe', 'start_probe', 'end_probe', 'chr', 'start', 'end', sample)]
  df <- df[!duplicated(df[,c('chr','start','end')]),]
  if(sample==sample_list[1]) {
    df_all <- df
  } else {
    df_all <- merge(df_all, df, by=c('probe_id','chr_probe', 'start_probe', 'end_probe', 'chr', 'start', 'end'), all=TRUE)
  }
}
df_all$Mean_meth_across_samples <- NA
df_all[nrow(df_all) + 1,] <- NA
df_all[nrow(df_all), 7] <- 'Mean_meth_across_CpG'
for(j in 8:(ncol(df_all)-1)) {
  mean_value <- mean(as.numeric(as.vector(unlist(df_all[1:(nrow(df_all)-1),j]))), na.rm=TRUE)
  mean_value <- format(round(mean_value, 3), nsmall = 3)
  df_all[nrow(df_all), j] <- mean_value
}
for(i in 1:nrow(df_all)) {
  mean_value <- mean(as.numeric(as.vector(unlist(df_all[i, 8:(ncol(df_all)-1)]))), na.rm=TRUE)
  mean_value <- format(round(mean_value, 3), nsmall = 3)
  df_all[i, ncol(df_all)] <- mean_value
}
write.csv(df_all, paste(work_dir, file_name_meth, sep = ''), row.names = FALSE)


### for coverage
for(sample in sample_list) {
  df <- read.table(paste(work_dir, sample, file_appendix, sep = ''), header = FALSE)
  df <- df[,c(1,2,3,4,5,6,8)]
  colnames(df) <- c('chr_probe', 'start_probe', 'end_probe', 'chr', 'start', 'end', sample)
  df$probe_id <- paste(df$chr_probe, df$start_probe, df$end_probe, sep='_')
  df <- df[,c('probe_id','chr_probe', 'start_probe', 'end_probe', 'chr', 'start', 'end', sample)]
  df <- df[!duplicated(df[,c('chr','start','end')]),]
  if(sample==sample_list[1]) {
    df_all <- df
  } else {
    df_all <- merge(df_all, df, by=c('probe_id','chr_probe', 'start_probe', 'end_probe', 'chr', 'start', 'end'), all=TRUE)
  }
}
df_all$Mean_Coverage_across_samples <- NA
df_all[nrow(df_all) + 1,] <- NA
df_all[nrow(df_all), 7] <- 'Mean_coverage_across_CpG'
for(j in 8:(ncol(df_all)-1)) {
  mean_value <- mean(as.numeric(as.vector(unlist(df_all[1:(nrow(df_all)-1),j]))), na.rm=TRUE)
  mean_value <- format(round(mean_value, 3), nsmall = 3)
  df_all[nrow(df_all), j] <- mean_value
}
for(i in 1:nrow(df_all)) {
  mean_value <- mean(as.numeric(as.vector(unlist(df_all[i, 8:(ncol(df_all)-1)]))), na.rm=TRUE)
  mean_value <- format(round(mean_value, 3), nsmall = 3)
  df_all[i, ncol(df_all)] <- mean_value
}
write.csv(df_all, paste(work_dir, file_name_coverage, sep = ''), row.names = FALSE)






