#Genome-Wide Association Study of Crohn’s Disease
Note: Files in this repository need to be downloaded to be viewed properly. The "Invalid PDF" message is just a GitHub display limitation.

library(data.table)




pheno <- read_table("pheno (1).txt")
head(pheno)
fam_data <- fread("final.fam")
colnames(fam_data) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")
bim_data <- fread("final.bim")
colnames(bim_data) <- c("Chromosome", "SNP", "GeneticDistance", "Position", "Allele1", "Allele2")
bed_data <- fread("final.bed")

#Quality Control
##missingness
###first iteration
plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",  
  "--bfile final",       
  "--geno 0.05",         # Keep SNPs with missingness <5%
  "--make-bed",          
  "--out final_geno05"   
)
system(plink_command)

plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",
  "--bfile final_geno05",  
  "--mind 0.10",           # Keep individuals with missingness <10%
  "--make-bed",            
  "--out final_mind10"     
)
system(plink_command)

###second iteration
plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",
  "--bfile final_mind10",  
  "--geno 0.01",           # Keep SNPs with missingness <1%
  "--make-bed",            
  "--out final_geno01"     
)
system(plink_command)

plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",
  "--bfile final_geno01",  
  "--mind 0.02",           # Keep individuals with missingness <2%
  "--make-bed",            
  "--out final_mind02"    
)
system(plink_command)

#perform missingness test
plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",
  "--bfile final_mind02",         
  "--test-missing",               # Perform missingness test
  "--out final_missing_test"      
)
system(plink_command)

final_missing <- fread("final_missing_test.missing")
head(final_missing)
filtered_snps <- final_missing[P < 0.05, .(SNP)]
head(filtered_snps)
# Save the SNPs to exclude into a text file
write.table(filtered_snps$SNP, "exclude_snps.txt", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)

plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",
  "--bfile final_mind02",         
  "--exclude exclude_snps.txt",   # Exclude SNPs based on the file created
  "--make-bed",                   
  "--out final_cleaned_no_missing"  
)
system(plink_command)

missing_clean <- fread("final_cleaned_no_missing.bim")

##MAF
plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",
  "--bfile final_cleaned_no_missing",  
  "--maf 0.01",                       # Filter SNPs with MAF >= 1%
  "--make-bed",                        
  "--out final_cleaned_no_missing_maf_01"  
)
system(plink_command)

##Hardy-Weinberg
plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",
  "--bfile final_cleaned_no_missing_maf_01",  
  "--hardy",                      
  "--make-bed",                        
  "--out final_cleaned_no_missing_maf_hwe"  
)
system(plink_command)
hard <- fread("final_cleaned_no_missing_maf_hwe.hwe")
hard2 = hard[hard$TEST=="UNAFF",]
summary(hard2$P)

#The p-values seem to be mostly above 0.05 (the usual threshold for statistical significance), 
#indicating that the SNPs generally do not significantly deviate from Hardy-Weinberg equilibrium 
#in the unaffected individuals.

hard2_all <- hard[hard$TEST == "ALL",]
summary(hard2_all$P)

# Filter SNPs where TEST is unaffected
hard_controls_filtered <- hard2[hard2$P > 1e-6, ]

# Filter SNPs where TEST is affected
hard_cases <- hard[hard$TEST == "AFF", ]
hard_cases_filtered <- hard_cases[hard_cases$P > 1e-10, ]

# Combine the filtered results
filtered_hwe <- rbind(hard_controls_filtered, hard_cases_filtered)
summary(filtered_hwe$P)

# Save the filtered HWE data to a file
write.table(filtered_hwe$SNP, "hwe_filtered_snps.txt", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
# Identify SNPs to exclude
excluded_snps <- setdiff(hard$SNP, filtered_hwe$SNP)

# Save the excluded SNPs
write.table(excluded_snps, "hwe_excluded_snps.txt", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",
  "--bfile final_cleaned_no_missing_maf_01",
  "--exclude hwe_excluded_snps.txt",  # Exclude failing SNPs
  "--make-bed",
  "--out final_cleaned_no_missing_maf_hwe_excluded"
)
system(plink_command)
#####
##Heterozygosity
plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",
  "--bfile final_cleaned_no_missing_maf_hwe_excluded",
  "--het",                       # Calculate heterozygosity
  "--out het_test"               
)
system(plink_command)
het <- fread("het_test.het")
summary(het$F)

# Visualize the F values
hist(het$F, main = "Distribution of Inbreeding Coefficients (F)", xlab = "F")
#the histogram shows that there seem to be no unsually low or high heterozygosity

#remove individuals with heterozygosity +/-3 SD from the mean
# Calculate mean and standard deviation
mean_F <- 0.001807  
sd_F <- sd(het$F, na.rm = TRUE)

# Define thresholds
lower_limit <- mean_F - 3 * sd_F  # Lower threshold
upper_limit <- mean_F + 3 * sd_F  # Upper threshold
# Identify outliers
outliers <- het$F < lower_limit | het$F > upper_limit

# Filter out outliers
filtered_het <- het[!outliers, ]

write.table(filtered_het, "filtered_heterozygosity_3SD.txt", quote = FALSE
            , row.names = FALSE)

######
##relatedness
plink_command <- paste("~/Documents/R/plink_mac_20241022/plink",
                       "--bfile final_cleaned_no_missing_maf_hwe_excluded",
                       "--genome",
  "--make-bed",
  "--out relatedness"
)
system(plink_command)

relatedness <- read.table("relatedness.genome", header = TRUE)
head(relatedness)

#identifying pairs with high-relatedness by filtering out pvalue> 0.5
high_relatedness <- relatedness[relatedness$PI_HAT > 0.5, ]
head(high_relatedness)
summary(high_relatedness$Z0)

#based on the summary stats, relatedness in the data set remains very low, meaning 
#that on average, individuals in the data set are not closely related

# Read bim file
bim_data <- fread("final_cleaned_no_missing_maf_hwe_excluded.bim")
nrow(bim_data)
# Check what chromosomes are present
table(bim_data$V1)  # V1 is the chromosome column

#no need for sex check

###Population Stratification Analysis

# Generate pairwise IBS matrix
plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",
  "--bfile final_cleaned_no_missing_maf_hwe_excluded",
  "--cluster",
  "--mds-plot 10",  # Calculate first 10 dimensions
  "--out population_struct"
)
system(plink_command)

pop_struct <- fread("population_struct.mds")
summary(pop_struct)

# Calculate variance explained by each dimension
mds_vars <- apply(pop_struct[,c("C1", "C2", "C3", "C4", "C5", "C6", "C7", 
                                "C8", "C9", "C10")], 2, var)
var_explained <- mds_vars/sum(mds_vars) * 100

# Plot variance explained
barplot(var_explained, 
        names.arg = paste0("C", 1:10),
        main = "Variance Explained by Each MDS Component",
        ylab = "Percent of Variance Explained")

# Plot first two dimensions
ggplot(pop_struct, aes(x = C1, y = C2)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Population Structure: First Two MDS Components",
       x = "Component 1",
       y = "Component 2")

# Plot second two dimensions
ggplot(pop_struct, aes(x = C2, y = C3)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Population Structure: First Two MDS Components",
       x = "Component 2",
       y = "Component 3")

##based on the above plots, we see a clear grouping of three clusters
##with some outliers towards the sides, component 1 and 2 capture 
##the most meaningful population structure compared to plot 2, so I will
##use them as covariates for the GWAS analysis


# Check sample sizes for each group
table(pheno$Disease)

# Run primary GWAS analysis with Crohn's Disease vs Control

# select the MDS components we want to use
pca_covariates <- pop_struct[, c("FID", "IID", "C1", "C2", "C3")]  # Using first 3 components

write.table(pca_covariates, 
            "pca_covariates.txt", 
            quote=FALSE, 
            row.names=FALSE)

# Create binary phenotype file for CD vs Control
pheno_cd <- pheno[pheno$Disease %in% c("CD", "Control"), ]
# Recode: Control = 1, CD = 2 
pheno_cd$Affectation <- ifelse(pheno_cd$Disease == "CD", 2, 1)

# Save phenotype file
write.table(pheno_cd[, c("Family_ID", "Indiv_ID", "Affectation")],
            "cd_pheno.txt", 
            row.names=FALSE, quote=FALSE)

# Run GWAS
## run a logit model since CD, the dependent variable, is binary,
## and also estimates probability of having a case
plink_command <- paste(
  "~/Documents/R/plink_mac_20241022/plink",
  "--bfile final_cleaned_no_missing_maf_hwe_excluded",
  "--pheno cd_pheno.txt",
  "--covar pca_covariates.txt",
  "--logistic",
  "--adjust",
  "--out cd_gwas_results"
)
system(plink_command)

gwas_results <- fread("cd_gwas_results.assoc.logistic")
head(gwas_results)

# Filter for only ADD tests
snp_results <- gwas_results[TEST=="ADD"]
significant_snps <- snp_results[P < 5e-8]
head(significant_snps)

# Visualize results
# Create Manhattan plot
#install.packages("qqman")
library(qqman)

manhattan_data <- gwas_results[TEST=="ADD", .(CHR=CHR, BP=BP, SNP=SNP, P=P)]

manhattan(manhattan_data,
          chr="CHR",
          bp="BP",
          p="P",
          snp="SNP",
          main="Manhattan Plot of CD GWAS",
          suggestiveline = -log10(1e-5),  # Less stringent threshold
          genomewideline = -log10(5e-8))  # Standard GWAS significance

## We see that at chr 1 and 16 there are significant peaks, implying strong signals
## at these points with the strongest at chr 16. We also see multiple peaks indicating the
## polygenic nature of Crohn's Disease. The multiple snps clustered together at these
## points indicate true genetic signals as this shows that the nearby SNPs are also correlated

# Create QQ plot to check for inflation
qq(manhattan_data$P, main="Q-Q plot of GWAS p-values")

## the qq plot further confirms that some SNPs do show true association
# Get significant SNPs (p < 5e-8)
top_snps <- gwas_results[TEST=="ADD" & P < 5e-8][order(P)]
# Extract SNP IDs
snp_ids <- top_snps$SNP
head(snp_ids)

# Load and check annotation file
annot <- fread("Ichip_anno_simplified_03042019 (2).txt")
head(annot$IlluminaID)

# Filter Chr 16 region of interest (~49.3Mb)
chr16_annot <- annot[chr == "16"]
chr16_region <- chr16_annot[grep("4931|4932|4930", IlluminaID)]

# Filter Chr 1 region of interest (~67.4Mb)
chr1_sig <- top_snps[CHR == 1]
chr1_annot <- annot[chr == "1"]
chr1_region <- chr1_annot[grep("6744|6740|6745", IlluminaID)]

# Get summary statistics for significant regions
# Chr 16
chr16_summary <- top_snps[CHR == 16, .(
  n_snps = .N,
  min_p = min(P),
  max_or = max(OR)
)]

# Chr 1
chr1_summary <- top_snps[CHR == 1, .(
  n_snps = .N,
  min_p = min(P),
  max_or = max(OR)
)]

#based on the summary statistics, we can see that there are 13 snps in chromosme 16
# and has an odds ratio of 1.81, which makes it a risk variant, increasing Crohn's Disease risk
# Chromosome 1 has 3 significant snps and has an odds ratio of 0.66


# Format  GWAS results for FUMA
fuma_format <- gwas_results[TEST=="ADD", .(
  SNP = SNP,
  CHR = CHR,
  BP = BP,
  A1 = A1,
  A2 = ifelse(exists("A2", where=gwas_results), A2, NA),  # Reference allele if available
  P = P,
  BETA = log(OR),  
  SE = STAT/log(OR), 
  N = NMISS
)]

head(fuma_format)
write.table(fuma_format, 
            "cd_gwas_for_fuma.txt", 
            row.names=FALSE, 
            quote=FALSE, 
            sep="\t")
