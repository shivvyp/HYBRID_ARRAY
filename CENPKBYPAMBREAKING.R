all_guides <-read.csv("/Users/shivanipatel/Desktop/Data/all_guide.csv")

library(dplyr)

# Add the "SNP count" column
filtered_grna <- all_guides %>% mutate(SNP_count = ifelse(lengths(strsplit(as.character(snp_posiions), ",")) == 1, "One", "Multi")) 
filtered_grna <- filtered_grna %>% mutate(pam_breaker = targets_ref != targets_alt)
                             
                             
                             
all_guides_ref = filtered_grna %>% dplyr::select(chrom=chrom, grna_start=grna_start, grna_end=grna_end,targets=targets_ref,pam_breaker=pam_breaker,strand=strand,snp_positions=snp_posiions,pam=pam_ref,
                                                 grna_seq_with_pam=grna_seq_ref,grna_seq=grna_seq_ref_t,notes=notes_ref,snps_in_seed=in_seed,lookup_ref=ref_lookup_ref,
                                                 lookup_alt=ref_lookup_alt) %>% mutate(intended_target="REF")
all_guides_alt = filtered_grna %>% dplyr::select(chrom=chrom, grna_start=grna_start, grna_end=grna_end,pam_breaker=pam_breaker,targets=targets_alt,strand=strand,snp_positions=snp_posiions,pam=pam_alt,
                                                 grna_seq_with_pam=grna_seq_alt,grna_seq=grna_seq_alt_t,notes=notes_alt,snps_in_seed=in_seed,lookup_ref=alt_lookup_ref,
                                                 
                                            lookup_alt=alt_lookup_alt) %>% mutate(intended_target="ALT")

all_guides_long <- rbind(all_guides_alt, all_guides_ref)
nrow(all_guides_long)
all_guides_long <- all_guides_long %>% filter(!grepl("PolyT", notes))
nrow(all_guides_long)
all_guides_long= all_guides_long %>% filter(!grepl("[ATCG]AG",pam))
nrow(all_guides_long)
all_guides_long

#filter look_ref and look_alt 

all_guides_long <-all_guides_long %>% mutate(total_sum = lookup_alt + lookup_ref)
all_guides_long <- all_guides_long %>% filter(total_sum == 1)
nrow(all_guides_long)
all_guides_long <- all_guides_long %>% filter(grepl(TRUE, targets))
nrow(all_guides_long)

#now subset the tables into PAM breakes
PAM_breaker_all <- all_guides_long %>% filter(grepl(TRUE, pam_breaker))
nrow(PAM_breaker_all)

#now take all the non-PAM breakers and break that into targeting alt and ref 
all_guides_nonpambreaker <- all_guides_long %>% filter(!pam_breaker & snps_in_seed)
nrow(all_guides_nonpambreaker)

#now take make a column which counts the number of snps within the guide, excluding the PAM and sort by the amount of gRNA's within the PAM  
# Define a function to count SNPs excluding positions 20, 21, and 23
count_snps <- function(snp_posiions) {
  # Split the SNP positions by comma and convert to a numeric vector
  snp_positions <- as.numeric(unlist(strsplit(snp_positions, ",")))
  
  # Exclude positions 20, 21, and 23 from the count
  snp_positions <- snp_positions[!(snp_positions %in% c(20, 21, 23))]
  
  # Count the remaining SNPs
  num_snps <- length(snp_positions)
  
  return(num_snps)
}



# Function to count SNPs in guide (excluding PAM)
count_snps <- function(snp_positions) {
  snp_positions <- as.numeric(unlist(strsplit(snp_positions, ",")))
  snp_positions <- snp_positions[!(snp_positions %in% c(22, 21, 23))]
  num_snps <- length(snp_positions)
  return(num_snps)
}


all_guides_nonpambreaker <- all_guides_nonpambreaker  %>%
  mutate(snp_count = sapply(snp_positions, count_snps))


#Function to count SNPs only in Seed region 
count_snps_in_seed <- function(snp_positions) {
  snp_positions <- as.numeric(unlist(strsplit(snp_positions, ",")))
  seed_start <- 13
  seed_end <- 20    
  snp_positions_in_seed <- snp_positions[snp_positions >= seed_start & snp_positions <= seed_end]
  num_snps_in_seed <- length(snp_positions_in_seed)
  return(num_snps_in_seed)
}


all_guides_nonpambreaker <- all_guides_nonpambreaker  %>%
  mutate(snp_count = sapply(snp_positions, count_snps))

all_guides_nonpambreaker <- all_guides_nonpambreaker %>% 
  mutate(snp_seed_count = sapply(snp_positions, count_snps_in_seed))

# Now filter for single snp in seed and remove all gRNA with a single SNP at position 13
nonpambreaking_sin13singlesnp <- all_guides_nonpambreaker %>% filter(snp_seed_count <= 1)
nonpambreaking_sin13singlesnp <- nonpambreaking_sin13singlesnp %>% filter(!grepl("13", snp_positions))
nonpambreaking_sin13singlesnp <- nonpambreaking_sin13singlesnp %>% arrange(desc(snp_count))

#filter guides to only have mutiple snps
nonpambreaking_sin13singlesnp_1 <- nonpambreaking_sin13singlesnp %>% filter(snp_count > 1)



#Now split into BY and CEN.PK table 
nonpambreaking_filtered_CEN <- nonpambreaking_sin13singlesnp_1 %>% filter(grepl("ALT", intended_target))
nrow(nonpambreaking_filtered_CEN) #618
nonpambreaking_filtered_BY <- nonpambreaking_sin13singlesnp_1 %>% filter(grepl("REF", intended_target))
nrow(nonpambreaking_filtered_BY)
#608

# now sample the single snp randomly to get to 6000
SINGLESNPINWHOLEGUIDE <- nonpambreaking_sin13singlesnp %>% filter(snp_count == 1)

#Now split into BY and CEN.PK table 
filtered_allelesp_CENPK <- SINGLESNPINWHOLEGUIDE%>% filter(grepl("ALT", intended_target))
filtered_allelesp_BY <- SINGLESNPINWHOLEGUIDE %>% filter(grepl("REF", intended_target))


#Now I want to filter my table to have 

paired_gRNAs_BY_to_CENPK <- data.frame(
  BY_Chromosome = character(0),
  BY_Start = numeric(0),
  BY_End = numeric(0),
  CENPK_Chromosome = character(0),
  CENPK_Start = numeric(0),
  CENPK_End = numeric(0),
  Distance = numeric(0),
  BasePairDifferences = numeric(0),
  grna_seq_with_pam_BY = character(0),  # Add the BY grna_seq_with_pam column
  grna_seq_with_pam_CENPK = character(0)  # Add the CENPK grna_seq_with_pam column
)

# Loop through each BY gRNA and find the closest CENPK gRNA on the same chromosome
for (i in 1:nrow(filtered_allelesp_BY)) {
  by_gRNA <- filtered_allelesp_BY[i, ]
  by_chromosome <- by_gRNA$chrom
  by_start <- by_gRNA$grna_start
  by_end <- by_gRNA$grna_end
  
  # Filter CEN.PK positions to match the same chromosome
  matching_cenpk_positions <- filtered_allelesp_CENPK[filtered_allelesp_CENPK$chrom == by_chromosome, ]
  
  if (nrow(matching_cenpk_positions) > 0) {
    # Calculate distances for all matching CEN.PK gRNAs on the same chromosome
    distances <- abs(matching_cenpk_positions$grna_start - by_start)
    closest_index <- which.min(distances)
    closest_gRNA <- matching_cenpk_positions[closest_index, ]
    distance <- distances[closest_index]
    
    # Calculate base pair differences
    bp_differences <- abs(by_start - closest_gRNA$grna_start)
    
    # Extract grna_seq_with_pam from both BY and CENPK
    grna_seq_with_pam_BY <- by_gRNA$grna_seq_with_pam
    grna_seq_with_pam_CENPK <- closest_gRNA$grna_seq_with_pam
    
    paired_gRNAs_BY_to_CENPK <- rbind(
      paired_gRNAs_BY_to_CENPK,
      data.frame(
        BY_Chromosome = by_chromosome,
        BY_Start = by_start,
        BY_End = by_end,
        CENPK_Chromosome = by_chromosome,
        CENPK_Start = closest_gRNA$grna_start,
        CENPK_End = closest_gRNA$grna_end,
        Distance = distance,
        BasePairDifferences = bp_differences,
        grna_seq_with_pam_BY = grna_seq_with_pam_BY,  # Populate the BY grna_seq_with_pam column
        grna_seq_with_pam_CENPK = grna_seq_with_pam_CENPK  # Populate the CENPK grna_seq_with_pam column
      )
    )
  }
}


ggplot(paired_gRNAs_BY_to_CENPK, aes(x = Distance)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(x = "Distance", y = "Frequency", title = "Histogram of Distances") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 1000, by = 100))











# Randomly sample pairs
set.seed(1)  # Set a seed for reproducibility
sampled_pairs <- paired_gRNAs_BY_to_CENPK %>% filter(Distance == 0) %>% sample_n(1018)

BYguidesingle <- sampled_pairs %>% select(grna_seq_with_pam_BY) %>% rename(grna_seq_with_pam_BY = "grna_seq_with_pam")
CENguidesingle <- sampled_pairs %>% select(grna_seq_with_pam_CENPK) %>% rename(grna_seq_with_pam_CENPK = "grna_seq_with_pam")

singles_long <- rbind(BYguidesingle,CENguidesingle)

# Load the dplyr package if you haven't already
library(dplyr)


print(sampled_df)

singles_long <- singles_long %>% mutate("type" = c("singleseed"))
nonpambreaking_sin13singlesnp_1 <- nonpambreaking_sin13singlesnp_1 %>% mutate("type" = c("multiseed"))
PAM_breaker_all <- PAM_breaker_all %>% mutate("type"= c("PAMbreaker"))


totalsample <- all_guides_long %>% right_join(singles_long, by= "grna_seq_with_pam") 
Tot_guides <- bind_rows(PAM_breaker_all, nonpambreaking_sin13singlesnp_1, totalsample)


write.csv(Tot_guides, file = "/Users/shivanipatel/Desktop/TOTALGRNASHIV.csv", row.names = FALSE)

