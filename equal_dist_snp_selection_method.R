######################################################################################################
#SELECT A NUMBER OF SNPs FROM A MAP FILE TO BE EQUALLY SPACED AND PROPORTIONALY TO CHROMOSOME LENGTH #
######################################################################################################
library(data.table)

#set the path to your working directory
dir <- file.path("your_path")

#Create a file that has these columns for your map file: chromosome, first_SNP_pos, last_SNP_pos, total_chr_distance (last_SNP_pos-first_SNP_pos),	SNPs/chrom (SNPs per chromosome)

#Read the file with the chromosome start and end position
chr_length<- read.table("chr_length.txt", header = TRUE)
head(chr_length)
#calculate (relative) total length of the chromosomes 
total_length<- sum(chr_length$total_chr_distance)
snp_in_LD_panel<- 300 #change this number to the required number of SNPs in your low-density (LD) panel  
#divide the desired number of SNPs you want to select (e.g. 300) in the LD panel with the total map length
index<- snp_in_LD_panel/total_length
#multiply the length of each chr with the index above and round to get the number of SNPs that are going to be selected from each chromosome according to its relative length
n_snp_per_chr<- round(index*chr_length$total_chr_distance, 0)
sum(n_snp_per_chr)    #check the total number of SNPs selected, should sum up to the desired number of SNPs in the LD panel  
n_snp_per_chr

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load the map file with the SNP IDs, chromosome number and physical position columns
snpmap <-read.table("C:/Users/s1899268/OneDrive - University of Edinburgh/Temp_Trout_CF_CK/rainbow_trout_snp_info.txt",header=T,sep="\t")
head(snpmap)

#split the map file above to chromosomes (the output from split() is a list)
chr_list = split(snpmap, snpmap$Chr)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Take each chromosome and divide it's length into equally distanced parts (equidistant positions), according to the number of snps we want to select (LD)
#These equidistant theoretical positions will be then used to find the nearest real position on the map  
b_chr<- list()
for (i in 1:length(chr_length$chr)) {
  b_chr[[i]]<- round(seq(chr_length$first_SNP[i], chr_length$last_SNP[i], length.out = n_snp_per_chr[i]), digits = 0)
}

#Check that the length of the list for each chromosome equals the number of SNPs we want to keep for this chromosome
for (i in 1:length(chr_length$chr)){
  print(length(b_chr[[i]]))
} 



#extract the BPPos column from the dataframe list 
#chr_list[["1"]][["BPPos"]], chr_list[["2"]][["BPPos"]], ... lapply applies a function to every list element. 
#In this case, the function is [, which takes rows as its first argument (we leave it blank to indicate all rows), columns as its second argument. 
#It returns the results in a list.
chr_bppos = lapply(chr_list, "[", , "BPPos")
chr_snpid = lapply(chr_list, "[", , "SNPID")

###***###***###***###***###***###***###***###***###***###***###***###***###***###***
#CREATE A TABLE WITH THE POSITIONS WE WANT TO KEEP FOR EACH CHROMOSOME

Merge_a_b<- list()
for (i in 1:length(chr_length$chr)) {
  a=data.table(Value=as.numeric(chr_bppos[[i]]))
  a[,merge:=Value]
  
  b=data.table(Value=b_chr[[i]])
  b[,merge:=Value]
  
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  
  Merge_a_b[[i]]=a[b,roll='nearest']
  print(length(unique(Merge_a_b[[i]][["Value"]]))) 
  print(length(unique(Merge_a_b[[i]][["merge"]])))
}

###***###***###***###***###***###***###***###***###***###***###***###***###***###***

#Find the nearest row of the theoritical positions to the real positions on the map file
#returns the row of the nearest neighbor position instead of the position
library(FNN)

#with this method instead of removing the selected row from the list (chr_bppos) I replace it with a 0 to keep the row number correct for each chr
#I end up having unique rows with exactly the number of SNPs I want to select
s1<- c()
for (i in 1:length(b_chr[[1]])) {
  s1[i]<- knnx.index(chr_bppos[["1"]], b_chr[[1]][i], k=1)
  chr_bppos[["1"]][s1[i]] <- 0
  s1<- sort(s1)
  print(length(unique(s1)))
}
s2<- c()
for (i in 1:length(b_chr[[2]])) {
  s2[i]<- knnx.index(chr_bppos[["2"]], b_chr[[2]][i], k=1)
  chr_bppos[["2"]][s2[i]] <- 0
  s2<- sort(s2)
  length(unique(s2))
}
s3<- c()
for (i in 1:length(b_chr[[3]])) {
  s3[i]<- knnx.index(chr_bppos[["3"]], b_chr[[3]][i], k=1)
  chr_bppos[["3"]][s3[i]] <- 0
  s3<- sort(s3)
  length(unique(s3))
}
s4<- c()
for (i in 1:length(b_chr[[4]])) {
  s4[i]<- knnx.index(chr_bppos[["4"]], b_chr[[4]][i], k=1)
  chr_bppos[["4"]][s4[i]] <- 0
  s4<- sort(s4)
  length(unique(s4))
}
s5<- c()
for (i in 1:length(b_chr[[5]])) {
  s5[i]<- knnx.index(chr_bppos[["5"]], b_chr[[5]][i], k=1)
  chr_bppos[["5"]][s5[i]] <- 0
  s5<- sort(s5)
  print(length(unique(s5)))
}
s6<- c()
for (i in 1:length(b_chr[[6]])) {
  s6[i]<- knnx.index(chr_bppos[["6"]], b_chr[[6]][i], k=1)
  chr_bppos[["6"]][s6[i]] <- 0
  s6<- sort(s6)
  print(length(unique(s6)))
}
s7<- c()
for (i in 1:length(b_chr[[7]])) {
  s7[i]<- knnx.index(chr_bppos[["7"]], b_chr[[7]][i], k=1)
  chr_bppos[["7"]][s7[i]] <- 0
  s7<- sort(s7)
  print(length(unique(s7)))
}
s8<- c()
for (i in 1:length(b_chr[[8]])) {
  s8[i]<- knnx.index(chr_bppos[["8"]], b_chr[[8]][i], k=1)
  chr_bppos[["8"]][s8[i]] <- 0
  s8<- sort(s8)
  print(length(unique(s8)))
}
s9<- c()
for (i in 1:length(b_chr[[9]])) {
  s9[i]<- knnx.index(chr_bppos[["9"]], b_chr[[9]][i], k=1)
  chr_bppos[["9"]][s9[i]] <- 0
  s9<- sort(s9)
  print(length(unique(s9)))
}
s10<- c()
for (i in 1:length(b_chr[[10]])) {
  s10[i]<- knnx.index(chr_bppos[["10"]], b_chr[[10]][i], k=1)
  chr_bppos[["10"]][s10[i]] <- 0
  s10<- sort(s10)
  print(length(unique(s10)))
}



#create a list with the vectors of snps to keep for each chromosome
s<- list(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, s29)
#check the number of snps selected
k<- c()
for (i in 1:length(chr_length$chr)){
  k[i]<-length(unique(s[[i]]))
}
# 300 snps
print(sum(k))

#Unlist the SNP position list of the SNPs chosen for the LD panel
LD_pos<- unlist(b_chr)
#mysnps_ld<-snpmap[snpmap$BPPos %in%LD_pos,]
#Create a list with the SNP IDs to keep within each chromosome 
LD_snp_ids<- list()
for (i in 1:length(chr_length$chr)){
  LD_snp_ids[[i]]<- chr_snpid[[i]][s[[i]]]
}
LD_ids<- unlist(LD_snp_ids)

write.table(LD_ids, file="300_snps.txt" , quote=FALSE, col.names = FALSE, row.names = FALSE)
