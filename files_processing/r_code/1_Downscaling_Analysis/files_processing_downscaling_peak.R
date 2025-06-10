library(tidyverse)

setwd('~/')


### This table has an unique identifier for each file analyzed. 
table = read.table('Table_downsize_factor_ENCODE_bysample_threshold_250_0.9.txt', sep = ',', header = T )


setwd('~/counts_dir')


## Peak calling function. ---- 
##### It uses the count files for each experiment using the bed.file pressent in the Annotations folder

peak_calling_function = function(data, size){
  ### You just put the one data frame directly
  #### Default size = 152800. Filtering out transcripts larger than this size. 
  
  size1 = ifelse(size != 152800, size, 152800)
  
  
  non_snar = data.frame(Index2 = c(78959,78960,79280,79283,79286,79288,79290,79292,79294, 79296,79298))
  
  
  data = data[!data$V8 %in% non_snar$Index2,]
  data$V11 = data$V3 -  data$V2
  
  ### Index = V8
  ### Bin = V9
  ### Counts = V10
  data$V8 = as.factor(data$V8)
  
  data1 = data %>% group_by(V8) 
  data1 = data1[data1$V9 == 'Normal' & data1$V11 < size1,]
  data1$V8 = as.factor(data1$V8)
  data1$V10 = as.integer(data1$V10)
  data2 = data[data$V8 %in% data1$V8,]
  
  
  rm(data1, data)
  
  data2 = data2 %>% select(V1, V4:V11) 
  data2 = data2 %>% group_by(V1,V4,V5,V6,V7,V8,V9) %>% summarise_if(is.integer, sum)

  ## Excluding nearby regions that do not have the desired size
  data2[data2$V9 == 'K' & data2$V11 == 500,9] = NA
  data2[data2$V9 == 'KK' & data2$V11 == 5000,9] = NA
  data2[data2$V9 == 'KKK' & data2$V11 == 50000,9] = NA
  
  ## Normal represents the counts for the gene of interest
  datalength = data2[data2$V9 == 'Normal',c(6,7,9)]
  
  data_pois = data2
  
  rm(data2)
  #### Dividing only for the length of the bin
  
  ### Whole genome background
  frac_total = sum(data_pois[data_pois$V9 == 'Normal', 8], na.rm = T) /2913022398
  
  ### Estimating Lambdas
  data_pois$V12 = ifelse(data_pois$V9 == 'K' & data_pois$V11 == 1000, data_pois$V10/1000,
                         ifelse(data_pois$V9 == 'KK' & data_pois$V11 == 10000, data_pois$V10/10000,
                                ifelse(data_pois$V9 == 'KKK' & data_pois$V11 == 100000, data_pois$V10/100000, data_pois$V10 )))
  
  data_pois1 = data_pois
  data_pois1 = data_pois1 %>% select(-V11, -V10) 
  data_pois1 = data_pois1 %>% pivot_wider(names_from = V9, values_from = V12) 
  data_pois1 = data_pois1 %>% arrange(V8)  
  rm(data_pois)
  
data_pois1$LNormal = frac_total
  
  data_pois1 = data_pois1 %>% select(1:6, 10, 7,8,9,11)
  
  data_pois1 = data_pois1 %>% left_join(datalength, by = 'V8')
  
  data_pois1[,8:11] = apply(data_pois1[,8:11], 2, function(x) x*data_pois1$V11)
  
  #### Selecting Maximum Lambda for each window of our gene of interest
  data_pois1$MaxL = apply(data_pois1[,8:11], 1, max, na.rm = T)
  
  #### Calculating pvalues using poisson distribution
  data_pois1 = data_pois1 %>% mutate(pvalue = ppois(lambda = MaxL, q = Normal, lower.tail = F, log.p = F))
  
  #### P-value adjustment
  data_pois1$padj = p.adjust(data_pois1$pvalue, method = 'fdr', n = length(data_pois1$pvalue))
  
  data_pois1 = data_pois1%>% ungroup() %>% select(V1,V4,V8, padj) %>% unite('V1', V1,V4,V8, sep = '&&')
  
  return(data_pois1)
  
}

## All files processing.

all = list.files(pattern = 'ureter')

pvalue_file = read.delim(all[1],header = F, sep = '\t' )
pvalue_file$V2 = as.integer(pvalue_file$V2)
pvalue_file$V3 = as.integer(pvalue_file$V3)
colnames(pvalue_file)[10] = 'V10'
pvalue_file$V10 = as.integer(pvalue_file$V10)


pval_test = peak_calling_function(data = pvalue_file, size = 152800)
pval_test = pval_test[,1]

non_snar = data.frame(Index2 = c(78959,78960,79280,79283,79286,79288,79290,79292,79294, 79296,79298))
pvalue_file = pvalue_file[!pvalue_file$V8 %in% non_snar$Index2,]
pvalue_file = pvalue_file[pvalue_file$V9 == 'Normal',]
pvalue_file = pvalue_file[,c(1,4,8,10)]
pvalue_file = pvalue_file %>% unite('V1', V1,V4,V8, sep = '&&')
compilation_counts = pvalue_file %>% select(V1)

## In this step we are downscaling reads for each unique experiment. 
## If the experiment has technical replicates counts for each window will be added.
## After counts have been added, downscaling by the respective factor will be conducted. 
## Once downscaling to the file has been performed, peak calling using the above function is performed. 

for (i in 1:nrow(table)) {
  files = list.files(pattern = table[i, 1])
  
  if (length(files) == 1) {
    print('There is one')
    file1 = read.delim(files[1], header = F, sep = '\t')
    file1 = file1 %>% mutate(V10 = floor(file1$V10 * table[i, 5]))
    
    
    write.table(
      x = file1,
      paste0('Downgraded_', table[i, 1],  '_', table[i,3],'.bed', sep = ''),
      row.names = F,
      sep = '\t',
      quote = F
    )
  
   xfile = file1
   
   xfile$V2 = as.integer(xfile$V2)
   xfile$V3 = as.integer(xfile$V3)
   colnames(xfile)[10] = 'V10'
   xfile$V10 = as.integer(xfile$V10)
   
   xfile2 = peak_calling_function(xfile, 152800)
   colnames(xfile2)[2] = paste0(table[i,1], '_', table[i,3])
   
   non_snar = data.frame(Index2 = c(78959,78960,79280,79283,79286,79288,79290,79292,79294, 79296,79298))
   xfile = xfile[!xfile$V8 %in% non_snar$Index2,]
   xfile = xfile[xfile$V9 == 'Normal',]
   xfile = xfile[,c(1,4,8,10)]
   xfile = xfile %>% unite('V1', V1,V4,V8, sep = '&&')
   colnames(xfile)[2] = paste0(table[i,1], '_', table[i,3])
   xfile[,2] = sapply(xfile[,2], as.numeric)
   
   compilation_counts = left_join(compilation_counts, xfile, by = 'V1')
   pval_test = left_join(pval_test, xfile2, by = 'V1')
      
  }
  else{
    collect = read.delim(files[1], header = F, sep = '\t')
    collect$V2 = as.character(collect$V2)
    collect$V3 = as.character(collect$V3)
    collect$V8 = as.character(collect$V8)
    collect = collect %>% group_by(V1, V2, V3, V4, V5, V6, V7, V8, V9) %>% select(-V10)
    print('There is more')
    for (j in 1:length(files)) {
      file1 = read.delim(files[j], header = F, sep = '\t')
      file1$V2 = as.character(file1$V2)
      file1$V3 = as.character(file1$V3)
      file1$V8 = as.character(file1$V8)
      file1 = file1 %>% group_by(V1, V2, V3, V4, V5, V6, V7, V8, V9)
      colnames(file1)[10] = paste0('col', j, sep = '')
      collect = cbind(collect,file1[,10])
    }
    collect = collect %>% ungroup()
    collect$final = rowSums(collect[,10:ncol(collect)])
    collect$final = as.numeric(collect$final)
    collect = collect %>% select(V1, V2, V3, V4, V5, V6, V7, V8, V9, final)
    collect = collect %>% mutate(final = floor(collect$final * table[i,5]))
    
    write.table(
      x = collect,
      paste0('Downgraded_', table[i, 1], '_', table[i,3], '.bed', sep = ''),
      row.names = F,
      sep = '\t',
      quote = F
    )
    
    xfile = collect
    
    xfile$V2 = as.integer(xfile$V2)
    xfile$V3 = as.integer(xfile$V3)
    colnames(xfile)[10] = 'V10'
    xfile$V10 = as.integer(xfile$V10)
    
    xfile2 = peak_calling_function(xfile, 152800)
    colnames(xfile2)[2] = paste0(table[i,1], '_', table[i,3])
    
    non_snar = data.frame(Index2 = c(78959,78960,79280,79283,79286,79288,79290,79292,79294, 79296,79298))
    xfile = xfile[!xfile$V8 %in% non_snar$Index2,]
    xfile = xfile[xfile$V9 == 'Normal',]
    xfile = xfile[,c(1,4,8,10)]
    xfile = xfile %>% unite('V1', V1,V4,V8, sep = '&&')
    colnames(xfile)[2] = paste0(table[i,1], '_', table[i,3])
    xfile[,2] = sapply(xfile[,2], as.numeric)
    
    compilation_counts = left_join(compilation_counts, xfile, by = 'V1')
    pval_test = left_join(pval_test, xfile2, by = 'V1')
    
  }
}



write.table(pval_test, file = 'ATAC_ENCODE_TISSUE_Compilation_pValue_Gene_Downgraded_FINAL.txt', row.names = F, sep = ',', quote = F)

write.table(compilation_counts, file = 'ATAC_ENCODE_TISSUE_rawcounts_compilation_Downgraded_FINAL.txt', row.names = F, sep = ',', quote = F)
