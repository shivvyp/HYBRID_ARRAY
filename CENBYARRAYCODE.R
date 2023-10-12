
library(dplyr)

# Merge the data frames as shown in the previous response
merged_df <- data.frame(
  Merged_Column = paste(Tot_guides$chrom, Tot_guides$grna_start, Tot_guides$grna_end, sep = "_"),
  stringsAsFactors = FALSE
)

# Extract the "grna" column from each data frame and create a new data frame
grna_df <- data.frame(
  grna = c(Tot_guides$grna_seq),
  stringsAsFactors = FALSE
)

# Combine the "Name" and "grna" columns into a new data frame
final_df <- data.frame(
  Name = merged_df$Merged_Column,
  Sequence = grna_df$grna,
  stringsAsFactors = FALSE
)




####
library(Biostrings)
complete_set = final_df
#
pbs = "gcagtgaaagataggtgacc"
structural_region_with_binding = "gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtggtgc"
rev_primer="cactatagggcgaattgggtaccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGACATAAAAAACAAAAAAAgcaccaccgactcgg"
nchar(rev_primer)
oligos = paste(pbs,complete_set$Sequence,structural_region_with_binding,sep="")
final.oligoD= DNAStringSet(oligos)
Acnt=letterFrequency(final.oligoD, 'A')
Tcnt=letterFrequency(final.oligoD, 'T')
final.oligoDrc=final.oligoD
for(i in 1:length(final.oligoD) ){
  # print(i)
  if(Acnt[i]>Tcnt[i]){ final.oligoDrc[i]=reverseComplement(final.oligoDrc[i]) }
}
data.frame(Name=complete_set$Name,Sequence=final.oligoDrc) %>% write.csv("/Users/shivanipatel/Desktop/Total gRNA counts /oct1023_S_twist.csv")
#data.frame(name=complete_set$Name,sequence=final.oligoDrc,complete_set) %>% write.csv("oct1023_S_twist.csv")
### Check for BSTEII ####


#read in the file gain 
gRNAtable_full <- read.csv("/Users/shivanipatel/Desktop/Total gRNA counts /oct1023_S_twist.csv")

gRNAtable_full <- gRNAtable_full %>% rename(Sequence = "grna_seq_with_pam")

grna <- bind_cols(gRNAtable_full, Tot_guides, by = "grna_seq_with_pam")

write.csv(grna, file = "/Users/shivanipatel/Desktop/Total gRNA counts /wide_oct1023_S_.csv", row.names = FALSE)
