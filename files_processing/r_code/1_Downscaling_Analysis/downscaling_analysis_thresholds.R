#Packages ----
library(tidyverse



#Directories -----

#### if at home simonlizarazo
setwd('/Users/simonlizarazo/Dropbox/UIUC/1_VanBortle Lab/1_Enter_Here/')
#### if at lab vanbortlelab
#### all this files are available upon request

#### These files include the total number of genes considered as significant at each downscaling factor after applying the peak_calling function

setwd('/Users/vanbortlelab/Dropbox/UIUC/1_VanBortle Lab/1_Enter_Here/')


encode = read.csv('downsample_file_encode.csv')
encode = encode[,c(-1)]
encode$dataframe = 'ENCODE Data'
colnames(encode)[3] = 'sample'

brain = read.csv('downsample_file_brain.csv')
brain = brain |> unite('sample', sample, id, sep = '_')
brain = brain[,c(-1)]
brain$dataframe = 'Brain Data'

brain1 = read.csv('downsample_file_brain_biologicalrep.csv')
brain1$sample = gsub('_\\d+$', '', brain1$sample)
brain1$sample = gsub("_\\d+(\\.\\d+)?$", '', brain1$sample)

brain1 = brain1 |> unite('sample', sample, id, sep = '_')
brain1 = brain1[,-1]
brain1$dataframe = 'Brain Data (Biological Replicates)'


tcga = read.csv('downsample_file_TCGA.csv')
tcga = tcga[,-1]
colnames(tcga)[3] = 'sample'
tcga$dataframe = 'TCGA Data'

liver = read.csv('downsample_file_liver.csv')
liver = liver[,c(-1,-6)]
liver$dataframe = 'Liver Data'

all = bind_rows(encode, brain, tcga, liver, brain1)

all = all |> filter(Seqdepth > 0)

write.table(all, 'downsampling_analysis_all_datasets_genes_vs_seqdepth.txt', sep = ',')

### regression analysis
######## Let's fit a function
### Michaelis mendel function can fit this case (ymax*x/(x +a))
#### I will fit the data using sequence depth up to 2e9

#### Considering the lineweaver burk linearization
aa = all
aa = aa |> transmute(RNA.C = 1/RNA.Central, RNA.P = 1/RNA.POL.III, propR = 1 /prop_RNA, propP = 1/prop_pol3, Seqdepth = 1/Seqdepth)

plot6 = all |> ggplot() + geom_line(
  aes(y = RNA.Central, x = Seqdepth, color = sample),
  show.legend = F,
  alpha = 0.5
) + theme_minimal() + ggtitle('Global Representation') +
  theme(
    panel.border = element_rect(colour = 'black', fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(),
    aspect.ratio = 1
  ) + ylab('Number of significant genes') + xlab('Sequencing Depth') +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                     limits =  c(0 , 1e+09)) + geom_vline(aes(xintercept = 250000000),
                                                          color = 'black',
                                                          linetype = 'dashed')

y### Only using the real data, without the Brain data grouped by biological replciate
data1 = aa[-7603:-7911,]

a1 = lm(RNA.C~Seqdepth, data = data1)
a2 = lm(RNA.P~Seqdepth, data = data1)
a3 = lm(propR~Seqdepth, data = data1)
a4 = lm(propP~Seqdepth, data = data1)

coeff = list(RNA.C = a1$coefficients,
             RNA.P = a2$coefficients,
             propR = a3$coefficients,
             propP = a4$coefficients)

coeff = as.data.frame(coeff)
coeff[nrow(coeff) + 1 ,] = 1/coeff[1,]
coeff[nrow(coeff) + 1 ,] = coeff[2,] * coeff[3,]

### Removing the biological replicates for the brain data for regression analysis
frequence_sample2 = all[-7603:-7911,]


model_RNA.C = nls(RNA.Central ~ MAX*Seqdepth/(K + Seqdepth), data = frequence_sample2,
                  start = list(MAX = coeff[3,1],
                               K = coeff[4,1]))


model_Pol3 = nls(RNA.POL.III ~ MAX*Seqdepth/(K + Seqdepth), data = frequence_sample2,
                 start = list(MAX = coeff[3,2],
                              K = coeff[4,2]))


model_propRNA = nls(prop_RNA ~ MAX*Seqdepth/(K + Seqdepth), data = frequence_sample2,
                    start = list(MAX = coeff[3,3],
                                 K = coeff[4,3]))

model_propPol3 = nls(prop_pol3 ~ MAX*Seqdepth/(K + Seqdepth), data = frequence_sample2,
                     start = list(MAX = coeff[3,4],
                                  K = coeff[4,4]))


coeff2 = as.data.frame(list(RNA.C = coef(model_RNA.C),
                            RNA.P = coef(model_Pol3),
                            propR = coef(model_propRNA),
                            propP = coef(model_propPol3)))

coeff2[nrow(coeff2) + 1,] =  coeff2[1,]*0.9
coeff2[nrow(coeff2) + 1,] = coeff2[3,]*coeff2[2,]/(coeff2[1,]-coeff2[3,])

### ONLY FOR THE BRAIN DATA WE ADDED BIOLOGICAL REPLICATES WITH THE GOAL OF KEEPING THEM IN OUR ANALYSIS DUE TO LOW SEQUENCING DEPTH. 
### Bringin back the biological replicates for the saturation curve

frequence_sample2 = all
sat_count = data.frame(y = predict(model_propRNA, newdata = frequence_sample2), x = frequence_sample2$Seqdepth, t = frequence_sample2$dataframe)
sat_count = sat_count[order(sat_count$x),]
rownames(sat_count) = seq(1:nrow(sat_count))

sat_count1 = NULL
for(i in 1:nrow(sat_count)){
  yty = frequence_sample2 |> filter(Seqdepth >= sat_count[i,2])
  yty1 = yty |> filter(prop_RNA < 0.75)

  yty = yty |> filter(!sample %in% yty1$sample)

  sss = data.frame(count = sat_count[i,2], saturation = sat_count[i,1], samples = length(unique(yty[yty$dataframe != 'Brain Data (Biological Replicates)',3])),
                   encode = length(unique(yty[yty$dataframe == 'ENCODE Data',3])),
                   brain = length(unique(yty[yty$dataframe == 'Brain Data',3])),
                   brain_bio = length(unique(yty[yty$dataframe == 'Brain Data (Biological Replicates)',3])),
                   tcga = length(unique(yty[yty$dataframe == 'TCGA Data',3])),
                   liver = length(unique(yty[yty$dataframe == 'Liver Data',3])),
                   percentage = '0.75')


  sat_count1 = bind_rows(sat_count1, sss)

  yty = frequence_sample2 |> filter(Seqdepth >= sat_count[i,2])
  yty1 = yty |> filter(prop_RNA < 0.875)

  yty = yty |> filter(!sample %in% yty1$sample)



  sss = data.frame(count = sat_count[i,2], saturation = sat_count[i,1], samples = length(unique(yty[yty$dataframe != 'Brain Data (Biological Replicates)',3])),
                   encode = length(unique(yty[yty$dataframe == 'ENCODE Data',3])),
                   brain = length(unique(yty[yty$dataframe == 'Brain Data',3])),
                   brain_bio = length(unique(yty[yty$dataframe == 'Brain Data (Biological Replicates)',3])),
                   tcga = length(unique(yty[yty$dataframe == 'TCGA Data',3])),
                   liver = length(unique(yty[yty$dataframe == 'Liver Data',3])),
                   percentage = '0.875')

  sat_count1 = bind_rows(sat_count1, sss)

  yty = frequence_sample2 |> filter(Seqdepth >= sat_count[i,2])
  yty1 = yty |> filter(prop_RNA < 0.90)

  yty = yty |> filter(!sample %in% yty1$sample)

  sss = data.frame(count = sat_count[i,2], saturation = sat_count[i,1], samples = length(unique(yty[yty$dataframe != 'Brain Data (Biological Replicates)',3])),
                   encode = length(unique(yty[yty$dataframe == 'ENCODE Data',3])),
                   brain = length(unique(yty[yty$dataframe == 'Brain Data',3])),
                   brain_bio = length(unique(yty[yty$dataframe == 'Brain Data (Biological Replicates)',3])),
                   tcga = length(unique(yty[yty$dataframe == 'TCGA Data',3])),
                   liver = length(unique(yty[yty$dataframe == 'Liver Data',3])),
                   percentage = '0.90')


  sat_count1 = bind_rows(sat_count1, sss)
}
sat_count1 = sat_count1 |> distinct()
frequence_sample2 = frequence_sample2 |> distinct(Seqdepth, .keep_all = T)

frequence_sample2 = frequence_sample2 |> right_join(sat_count1[,c(1,3:9)], join_by('Seqdepth' =='count')) |>  distinct()


strips = strip_themed(background_x = elem_list_rect(fill =c('#ebbab9','#bf5e5c', '#910b09')),
                      text_x = elem_list_text(colour = 'black', face = 'bold'))

### Only regression line and Global
### All the samples
plot6.1 = frequence_sample2 |> filter(sample != 'Brain Data (Biological Replicates)') |> ggplot(aes(x = Seqdepth)) +
  geom_line(aes(y = prop_RNA, color = sample), show.legend = F, alpha = 0.01) +
  geom_line(aes(y = predict(model_propRNA, newdata = frequence_sample2), color = 'Predictive Model'), linewidth = 1) +
  geom_line( aes(y = samples/699, color = percentage), linewidth = 1 , show.legend = T) +
  facet_wrap2(~percentage, strip = strips, axes = 'all') + xlim(0,1e+09) +
  xlim(0,1e+09) +
  theme_bw() +
  scale_y_continuous( name = 'Significant Proportion of Genes',
                                   sec.axis = sec_axis(trans =  ~.*699, name = 'Samples Above Threshold')) +
  scale_color_manual(values = c('Predictive Model' ='black', '0.75'= '#ebbab9', '0.875' = '#bf5e5c', '0.90' = '#910b09')) +
  theme(text = element_text(family = 'Arial'),
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "#9e1818"),legend.title = element_blank()) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()), limits =  c(0 ,1e+09))


plot7 = frequence_sample2 |> ggplot(aes(x = Seqdepth)) +
  geom_line( aes(y = samples/699, color = 'Global'), linewidth = 1 ) +
  geom_line( aes(y = encode/156, color = 'ENCODE Data'), linewidth = 1 ) +
  geom_line( aes(y = brain/115, color = 'Brain Data'), linewidth = 1 ) +
  geom_line( aes(y = tcga/408, color = 'TCGA Data'), linewidth = 1 ) +
  geom_line( aes(y = liver/20, color = 'Liver Data'), linewidth = 1 ) +
  geom_line( aes(y = brain_bio/28, color = 'Brain Data (Biological Replicates)'), linewidth = 1 ) +
  xlim(0,1e+09) +
  theme_bw() +
  facet_wrap2(~percentage, strip = strips, axes = 'all') +
  theme(
    axis.title.y = element_text(color = "black")) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()), limits =  c(0 ,1e+09)) +
  scale_color_manual(values = c('Global'= '#9e1818', 'ENCODE Data' = '#f79525', 'Brain Data' = '#12b335',
                                'Brain Data (Biological Replicates)' = '#c5ebcf',
                                'TCGA Data' = '#0e7eb3', 'Liver Data' = '#3310b3' ),
                     breaks = c('Predictive Model', 'Global', 'TCGA Data','ENCODE Data', 'Brain Data', 'Brain Data (Biological Replicates)', 'Liver Data'))



### Samples by source
freq3 = frequence_sample2 |> pivot_longer(cols = samples:liver)
freq3$name = gsub('samples', 'Global', freq3$name)
freq3$name = gsub('tcga', 'TCGA', freq3$name)
freq3$name = gsub('encode', 'ENCODE', freq3$name)
freq3$name = gsub('brain_bio', 'Brain (Grouping)', freq3$name)
freq3$name = gsub('brain', 'Brain (No grouping)', freq3$name)

freq3$name = gsub('liver', 'Liver', freq3$name)
freq3$name = factor(freq3$name, levels= c('ENCODE', 'TCGA', 'Liver', 'Brain (No grouping)', 'Brain (Grouping)', 'Global'))

plot8 = freq3 |> filter(percentage == '0.90') |> ggplot(aes(x = Seqdepth)) +
  geom_line(aes(y = value, color = name), show.legend = F) +
  facet_wrap( ~ name, scales = 'free', nrow = 2) + xlim(0, 1e+09) + theme_minimal() + geom_vline(aes(xintercept = 250000000),
                                                                                       color = 'black',
                                                                                       linetype = 'dashed') +
  theme(
    aspect.ratio = 1,
    panel.border = element_rect(colour = 'black', fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line()
  ) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                     limits =  c(0 , 1e+09)) + xlab('Sequencing Depth') +
  ylab('Number of experiments included') 


ggarrange(a = ggarrange(plot1,plot2,plot3,plot4, plot5, plot6, nrow = 3, ncol = 2),
          plot8, nrow = 2, labels = c('A','B'), heights = c(1,1))

#### Plot with Final point for optimal depth sequence

plot6.2 = frequence_sample2 |> filter(sample != 'Brain Data (Biological Replicates)') |> ggplot(aes(x = Seqdepth)) +
  geom_line(aes(y = prop_RNA, color = sample), show.legend = F, alpha = 0.01) +
  geom_line(aes(y = predict(model_propRNA, newdata = frequence_sample2), color = 'Predictive Model'), linewidth = 1) +
  geom_line( aes(y = samples/699, color = percentage), linewidth = 1 , show.legend = T) +
  geom_segment(aes(x = 250000000, xend = 250000000, y = 0, yend = 0.9663878), colour = 'gray30', linetype = 'dashed') +
  geom_segment(aes(x = 0, xend = 250000000, y = 0.9663878, yend = 0.9663878), colour = 'gray30', linetype = 'dashed') +
  geom_point(aes(x = 250000000, y = 0.9663878), colour = 'gray30', size =3) +
  facet_wrap2(~percentage, strip = strips, axes = 'all') + xlim(0,1e+09) +
  xlim(0,1e+09) +
  theme_classic() +
  scale_y_continuous( name = 'Significant Proportion of Genes',
                      sec.axis = sec_axis(trans =  ~.*699, name = 'Samples Above Threshold')) +
  scale_color_manual(values = c('Predictive Model' ='black', '0.75'= '#ebbab9', '0.875' = '#bf5e5c', '0.90' = '#910b09')) +
  theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "#9e1818"),legend.title = element_blank()) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()), limits =  c(0 ,1e+09))

plot6.2

write.table(frequence_sample2, 'downsampling_analysis_regresion_analysis.txt', sep = ',')

#### Making tables
frequence_sample2 = all
ppp1 = frequence_sample2 |> ungroup() |> filter(Seqdepth >= 250000000)
pp = ppp1 |> filter(prop_RNA < 0.9)

ppp1 = ppp1 |> filter(!sample %in% pp$sample)
ppp1 = ppp1 |> filter(factor == 1.0)|> mutate( limit = 250000000)
ppp1$factor = ppp1$limit/ppp1$Seqdepth
ppp1$proof = ppp1$factor *ppp1$Seqdepth
ppp1 = ppp1 |> dplyr::select(3,4,7,8,9,10)
ppp1$dataframe = gsub(' ', '_', ppp1$dataframe)
ppp1$sample = gsub('_\\d+$', '', ppp1$sample)


aa = unique(ppp1$dataframe)

pp = ppp1 |> filter(dataframe == aa[1])

# This files are the ones that were used in the analysis passing our thresholds.

pa <- read.csv("~/Dropbox/UIUC/1_VanBortle Lab/1_Enter_Here/Sequencing_Methods/ATAC-Seq/ENCODE/Table_downsize_factor_ENCODE_bysample_FINAL.txt")

pp = pp |> left_join(pa[,c(1,3)], join_by('sample' == 'Experiment'))

pp = pp |> dplyr::select(1,3,7,5,2,6, 4)

write.table(pp,'Sequencing_Methods/ATAC-Seq/Table_downsize_factor_ENCODE_bysample_threshold_250_0.9.txt', quote = F, sep = ',' )

###
pp = ppp1 |> filter(dataframe == aa[2])

pp = pp |> dplyr::select(1,3,5,2,6, 4)
pp$sample = gsub("\\.", "-", pp$sample)

write.table(pp,'Sequencing_Methods/ATAC-Seq/Table_downsize_factor_TCGA_bysample_threshold_250_0.9.txt', quote = F, sep = ',' )

###
pp = ppp1 |> filter(dataframe == aa[3])

pp = pp |> dplyr::select(1,3,5,2,6, 4)

write.table(pp,'Sequencing_Methods/ATAC-Seq/Table_downsize_factor_LIVER_bysample_threshold_250_0.9.txt', quote = F, sep = ',' )

###
pp = ppp1 |> filter(dataframe == aa[4])

pp = pp |> dplyr::select(1,3,5,2,6, 4)

write.table(pp,'Sequencing_Methods/ATAC-Seq/Table_downsize_factor_BRAIN_bysample_threshold_250_0.9.txt', quote = F, sep = ',' )