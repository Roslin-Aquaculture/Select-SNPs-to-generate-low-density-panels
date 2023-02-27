######################################################################################################
#SELECT A NUMBER OF SNPs FROM A MAP FILE TO BE EQUALLY SPACED AND PROPORTIONALY TO CHROMOSOME LENGTH #
######################################################################################################
library(data.table)
library(FNN)

#Set the path to your working directory
dir <- file.path("your_path")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create a file that has the following columns in your text file: chromosome, first_SNP_pos, last_SNP_pos, total_chr_distance (last_SNP_pos-first_SNP_pos), SNPs/chrom (SNPs per chromosome)
#Read the file with the chromosome start and end position
chr_length<- read.table("chr_info.txt", header = TRUE)
head(chr_length)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Calculate the available total length of the chromosomes 
total_length<- sum(chr_length$total_chr_distance)
snps_in_LD_panel<- 200 #change this number to the required number of SNPs in your low-density (LD) panel  
#Divide the desired number of SNPs you want to select (e.g. 300) in the LD panel with the total map length
index<- snps_in_LD_panel/total_length
#Multiply the length of each chr with the index above and round to get the number of SNPs that are going to be selected from each chromosome according to its relative length
n_snps_per_chr<- round(index*chr_length$total_chr_distance, 0)
sum(n_snps_per_chr)    #check the total number of SNPs selected, should sum up to the desired number of SNPs in the LD panel  
n_snps_per_chr    #prints the number of SNPs to be selected from each chromosome according to its length

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Load the map file with the SNP IDs, chromosome number and physical position columns
snpmap <-read.table("snp_info.txt",header=T,sep="\t")
head(snpmap)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Split the map file above to chromosomes
chr_list = split(snpmap, snpmap$Chr)


#Take each chromosome and divide it's length into equally distanced parts (equidistant positions), according to the number of SNPs we want to select (snp_in_LD_panel)
#These equidistant theoretical positions will be then used to find the nearest real position on the map  
theor_pos<- list()
for (i in 1:length(chr_length$Chr)) {
  theor_pos[[i]]<- round(seq(chr_length$first_SNP[i], chr_length$last_SNP[i], length.out = n_snps_per_chr[i]), digits = 0)
}

#Check that the length of the list for each chromosome equals the number of SNPs we want to keep for this chromosome
for (i in 1:length(chr_length$Chr)){
  print(length(theor_pos[[i]]))
} 


#Extract the BPPos column from the dataframe list 
#In this case, the function is [, which takes rows as its first argument (we leave it blank to indicate all rows), columns as its second argument 
#It returns the results in a list
chr_bppos = lapply(chr_list, "[", , "BPPos")
chr_snpid = lapply(chr_list, "[", , "SNPID")

#Find the nearest row of the theoritical positions to the real positions on the map file
#Create a table with the theoretical positions we want to keep for each chromosome
Merge_real_theor<- list()
for (i in 1:length(chr_length$Chr)) {
  a=data.table(Value=as.numeric(chr_bppos[[i]]))
  a[,merge:=Value]
  
  b=data.table(Value=theor_pos[[i]])
  b[,merge:=Value]
  
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  
  Merge_real_theor[[i]]=a[b,roll='nearest']
  print(length(unique(Merge_real_theor[[i]][["Value"]]))) 
  print(length(unique(Merge_real_theor[[i]][["merge"]])))
}

#Returns the row of the nearest neighbor position instead of the position and we end up having a list with the SNPs we want to keep for each chromosome
s<- list()
for (x in 1:length(chr_length$Chr)){
  selected <- c()
  for (i in 1:length(theor_pos[[x]])) {
    selected[i]<- knnx.index(chr_bppos[[x]], theor_pos[[x]][i], k=1)
    chr_bppos[[x]][selected[[i]]] <- 0
    #selected<- sort(selected)
    #print(length(unique(selected)))
    s[[x]] <- selected
  }
}


#Check the total number of SNPs selected
k<- c()
for (i in 1:length(chr_length$Chr)){
  k[i]<-length(unique(s[[i]]))
}
print(sum(k))


#Unlist the SNP position list of the SNPs chosen for the LD panel
LD_pos<- unlist(theor_pos)
#mysnps_ld<-snpmap[snpmap$BPPos %in%LD_pos,]
#Create a list with the SNP IDs to keep within each chromosome 
LD_snp_ids<- list()
for (i in 1:length(chr_length$Chr)){
  LD_snp_ids[[i]]<- chr_snpid[[i]][s[[i]]]
}

LD_ids<- unlist(LD_snp_ids)

write.table(LD_ids, file="200_snps.txt" , quote=FALSE, col.names = FALSE, row.names = FALSE)
