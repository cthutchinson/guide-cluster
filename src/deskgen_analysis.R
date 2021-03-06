#install packages required if they haven't been already
list.of.packages <- c("readr", "tidyverse","stringr","Rtsne","parallel","dbscan","factoextra","ggpubr","repr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#install bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()

#load the required packages
library(readr)
library(tidyverse)
library(stringr)
library(Rtsne)
library(parallel)
library(org.Hs.eg.db)
library(GO.db)
library(dbscan)
library(factoextra)
library(ggpubr)
library(repr)

#load the data
example_guide_data <- read_delim("../data/example_guide_data.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
crispr_data = example_guide_data
View(crispr_data)

#create new columns to average the two repitions into one result for both drug and no drug conditions
crispr_data = mutate(crispr_data,avg_D7=(norm_count_D7_Rep1+norm_count_D7_Rep2)/2)
crispr_data = mutate(crispr_data,avg_D14=(norm_count_D14_Rep1+norm_count_D14_Rep2)/2)
crispr_data = mutate(crispr_data,avg_D7_PLX=(norm_count_PLX7_Rep1+norm_count_PLX7_Rep2)/2)
crispr_data = mutate(crispr_data,avg_D14_PLX=(norm_count_PLX14_Rep1+norm_count_PLX14_Rep2)/2)
#clean up columns
crispr_data = dplyr::select(crispr_data, gene_name,spacer_id,spacer_seq,plasmid=norm_count_plasmid,avg_D7,avg_D14,avg_D7_PLX,avg_D14_PLX)

#calculate the activity of the guide (-log2(count_t2/count_t1)) for each condition
crispr_data = mutate(crispr_data, D7_activity= -log2(avg_D7/plasmid), D7_activity_PLX = -log2(avg_D7_PLX/plasmid), D14_activity = -log2(avg_D14/plasmid), D14_activity_PLX = -log2(avg_D14_PLX/plasmid))

#function to count occurences of a character in a string (https://techoverflow.net/2012/11/10/r-count-occurrences-of-character-in-string/)
countCharOccurrences <- function(char, s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}

#count freqencies of bases in guide
crispr_data = mutate(crispr_data,A_freq = countCharOccurrences("A",spacer_seq)/20)
crispr_data = mutate(crispr_data,A_count = countCharOccurrences("A",spacer_seq))
crispr_data = mutate(crispr_data,T_freq = countCharOccurrences("T",spacer_seq)/20)
crispr_data = mutate(crispr_data,C_freq = countCharOccurrences("C",spacer_seq)/20)
crispr_data = mutate(crispr_data,G_freq = countCharOccurrences("G",spacer_seq)/20)
#include GC content 
crispr_data = mutate(crispr_data,GC_content = G_freq+C_freq)


#melting temperature according to Howley's formula 
#Tm = 64.9 + 41 * (nG + nC - 16.4) / (nA + nT + nG + nC))
crispr_data = mutate(crispr_data,temp_m = 64.9 + 41*(crispr_data$G_freq*20 + crispr_data$C_freq*20 - 16.4)/(crispr_data$A_freq*20 + crispr_data$T_freq*20 + crispr_data$G_freq*20 + crispr_data$C_freq*20))

#export guide sequences for analysis in bowtie
write_csv(as.data.frame(crispr_data$spacer_seq),"crispr_seq.txt",col_names=FALSE)

#./bowtie -v 0 hg19 -r reads/crispr_seq.txt > chrom.txt
#import results from bowtie matching to human reference genome
chrom = read_delim("../data/chrom.txt","\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
chrom = dplyr::select(chrom, chr=X3)
crispr_data = cbind(crispr_data,chrom)
crispr_data = mutate(crispr_data,chrom=str_replace(chr,"chr",""))
crispr_data = dplyr::select(crispr_data,-chr)

#export the sequences and run them in Oligoprop (MATLAB) to find the
#number of hairpins which can be formed in the sequence and the weights
num_hairpins = read_csv("../data/num_hairpins.txt",col_names = FALSE)
num_hairpins = as.integer(unlist(strsplit(as.character(num_hairpins[1,1]),split=NULL)))
crispr_data$hairpins = num_hairpins


#detect if sequence contains repetitive bases
GGGG = lapply(crispr_data$spacer_seq, function(x) str_detect(x,"GGGG"))
crispr_data$GGGG = as.numeric(GGGG)

#detect if the last base in the seqence is a G
G_end = lapply(crispr_data$spacer_seq, function(x) strsplit(x,split=NULL)[[1]][20]=="G")
crispr_data$G_end = as.numeric(G_end)
#detect if position 16 is a C
C_16 = lapply(crispr_data$spacer_seq, function(x) strsplit(x,split=NULL)[[1]][16]=="C")
crispr_data$C_16 = as.numeric(C_16)


#create a function to look up the biological function
annotate_term = function(x) {
   select(GO.db,filter(select(org.Hs.eg.db, keys=x, keytype="SYMBOL",columns="GO"),ONTOLOGY=="BP")[1,2],"TERM")[1,2]
}
 
#parallelize for efficiency
cl = makeCluster(detectCores()-1,type="FORK")
annotations = parSapply(cl,1:nrow(crispr_data), function(x) tryCatch({
  annotate_term(as.character(crispr_data[x,2]))
},error = function(e) {
  NA
}))

#annotating the genes took a long time, lets save them for reproduction sake
write_csv(as_tibble(annotations),path="gene_BP.txt",col_names = FALSE)

#read the gene annotations in from the bowtie file
annotations = read_csv("../data/gene_BP.txt", col_names = FALSE)
annotations = as_tibble(annotations)
annotations = dplyr::select(annotations,gene_BP=X1)
crispr_data = cbind(annotations,crispr_data)

#import deltaG from oligoprop
deltaG = read_csv("../data/deltaG.txt",col_names = FALSE)
deltaG = as.double(unlist(strsplit(as.character(deltaG[1,1]),split=" ")))
crispr_data$deltaG = deltaG

#tsne analysis
tsne_data = dplyr::select(crispr_data,A_freq,temp_m,hairpins,deltaG,GGGG,G_end,C_16,GC_content,chrom)
tsne = Rtsne(tsne_data, dims = 2, perplexity=400,verbose=TRUE,max_iter=500,check_duplicates=FALSE)
tsne_results = tsne$Y
tsne_results =as_tibble(tsne_results)
ggplot(tsne_results) + geom_point(mapping=aes(x=tsne_results[[1]],y=tsne_results[[2]],color=crispr_data$D14_activity)) + scale_colour_gradient(low = "blue", high = "yellow")

#cluster the t-sne output
db = dbscan(tsne_results,eps=1.5,minPts = 200)
fviz_cluster(db, data = tsne_results, stand = FALSE,
             ellipse = FALSE, show.clust.cent = FALSE,
             geom = "point",palette = "jco", ggtheme = theme_classic())
crispr_data$clusters = db$cluster
