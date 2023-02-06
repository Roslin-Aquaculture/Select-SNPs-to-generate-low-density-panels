# Select equaly spaced SNPs proportionally to chromosome length

This R script selects SNPs from a given map file and outputs a list with the SNP IDs and SNP positions to create a low-density panel of a predetermined size. The SNPs selected are equally distributed across the genome and proportionally to chromosome length. Additionally, this script always selects the first and the last SNP of each chromosome.  

##0. Input files
Two files are needed for this script:
1. The first file (chr_info.txt) is a text file with the following columns: Chr, first_SNP_pos, last_SNP_pos, total_chr_distance (last_SNP_pos-first_SNP_pos),	SNPs/chrom (SNPs per chromosome)
2. The second file (snp_info.txt) contains three columns: SNPID, Chr and BPPos

##1. The first part of the script prints out the number of SNPs that are going to be selected from each chromosome according to its length once we define the number of SNPs in the low-density panel.

```
total_length<- sum(chr_length$total_chr_distance)
snp_in_LD_panel<- 300 #change this number to the required number of SNPs in your low-density (LD) panel  
#Divide the desired number of SNPs you want to select (e.g. 300) in the LD panel with the total map length
index<- snp_in_LD_panel/total_length
#Multiply the length of each chr with the index above and round to get the number of SNPs that are going to be selected from each chromosome according to its relative length
n_snp_per_chr<- round(index*chr_length$total_chr_distance, 0)
sum(n_snp_per_chr)    #check the total number of SNPs selected, should sum up to the desired number of SNPs in the LD panel  
n_snp_per_chr
```

##2. In the second part we take each chromosome and divide it's length into equally distanced parts (equidistant positions), according to the number of SNPs we want to select (snp_in_LD_panel). These equidistant theoretical positions will be then used to find the nearest, real position on the map.  

```
b_chr<- list()
for (i in 1:length(chr_length$Chr)) {
  b_chr[[i]]<- round(seq(chr_length$first_SNP[i], chr_length$last_SNP[i], length.out = n_snp_per_chr[i]), digits = 0)
}
```

##3. In this step we find the nearest real positions on the map file of the theoritical positions, which we generated in the previous step.
(returns the coresponding rows of the nearest neighbor positions instead of the positions)

```
s<- list()
for (x in 1:length(chr_length$Chr)){
  test <- c()
  for (i in 1:length(b_chr[[x]])) {
  test[i]<- knnx.index(chr_bppos[[x]], b_chr[[x]][i], k=1)
  chr_bppos[[x]][test[[i]]] <- 0
  #test<- sort(test)
  #print(length(unique(test)))
  s[[x]] <- test
  }
  }
```

##4. Finally create a list with the equally spaced SNP IDs to keep within each chromosome according to its length

```
LD_snp_ids<- list()
for (i in 1:length(chr_length$Chr)){
  LD_snp_ids[[i]]<- chr_snpid[[i]][s[[i]]]
}
```