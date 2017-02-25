suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(stringr))


option_list <- list(
    make_option(c("-n", "--ncrna"),
                help = "analyze ncRNA instead of protein-coding loci",
                action = "store_true",
                default = FALSE),
    make_option(c("-o", "--output"),
                help = "output results to specified file name (by default, results are printed to terminal)",
                type = "character",
                default = "STDOUT")

)

parser <- OptionParser(option_list = option_list, description = "MMP app")

arguments <- parse_args(parser, positional_arguments = 1)

opts <- arguments$options

## Variables

# Directory containing mmp_mutations.txt
mmp_dir <- "/home/sfrenk/Documents/Resources/mmp"

# Strains to be removed from analysis due to missing data
filtered_strains = c("VC10116", "VC10118", "VC20007", "VC20009", "VC20495", "VC20507", "VC20539", "VC20543","VC40186", "VC40543", "VC20228", "VC20496", "VC20499", "VC20530", "VC20545")


## READ IN DATA

# Read in mmp strain data (this file lists all detected mutations in each mmp strain)
print("Reading in mmp mutation data...")
mmp <- read.table(paste0(mmp_dir, "/mmp_mutations.txt"), sep = "\t", header = TRUE)

# Read in input strain list
input_strains <- read.table(arguments$args)$V1
n_strains_input <- length(input_strains)

# Select interesting mutations
if (!opts$ncrna) {
    print("Finding protein coding genes...")
    mmp <- mmp[!(mmp$effect %in% c("synonymous", "ncRNA")) & ((mmp$feature == "coding_exon") | (mmp$effect == "splicing")),]

} else {
    print("Finding ncRNA loci...")
    mmp <- mmp[mmp$effect %in% c("ncRNA", "piRNA", "miRNA", ''),]

}

# Keep important variables and remove duplicates and strains with missing data
mmp <- mmp[c("strain", "CGC")]
mmp <- mmp[!duplicated(mmp),]
mmp <- mmp[!(mmp$strain %in% filtered_strains),]
n_strains_total <- length(unique(mmp$strain))


## FIND GENES

## Find genes present in input strains
mmp_input_genes <- plyr::count(mmp[mmp$strain %in% input_strains,], vars = "CGC")
mmp_input_genes <- mmp_input_genes[!is.na(mmp_input_genes$CGC),]
colnames(mmp_input_genes) <- c("gene", "input_count")

# Count occurances of these genes in the entire MMP dataset
mmp_total_genes <- plyr::count(mmp[mmp$CGC %in% mmp_input_genes$gene,], vars = "CGC")
mmp_total_genes <- mmp_total_genes[!is.na(mmp_total_genes$CGC),]
colnames(mmp_total_genes) <- c("gene", "total_count")
mmp_table <- merge(mmp_input_genes, mmp_total_genes)


## CALCULATE GENE ENRICHMENT USING FISHER'S EXACT TEST

mmp_hyper <- function(x, m, n){
    # x = "number of occurances" (number of time a particular gene was mutated)
    # m = number of times gene is mutated in mmp
    # n = number of times gene is wild-type in mmp
    # n_strains_input = size of subpopulation
    test <-phyper(x, m, n, n_strains_input, lower.tail = F)
    return(test)
}

# Note: have to subtract 1 from input gene count to determine P(x >= n)
mmp_table$pval_raw <- mapply(mmp_hyper,
                        mmp_table$input_count-1,
                        mmp_table$total_count,
                        (n_strains_total-mmp_table$total_count))


mmp_table$pval_adj <- p.adjust(mmp_table$pval_raw, "fdr", nrow(mmp_table))
mmp_table <- mmp_table[order(mmp_table$pval_raw),]

## GET GENE DESCRIPTIONS USING BIOMART
mart <- useMart(biomart="ensembl", dataset="celegans_gene_ensembl")

results <- getBM(attributes = c("external_gene_name", "description"), 
                 filters = "external_gene_name", 
                 values = mmp_table$gene,
                 mart = mart)
colnames(results) <- c("gene", "description")

results_table <- merge(mmp_table, results, by = "gene", all.x = TRUE)
results_table <- results_table[!(duplicated("gene")),]
results_table <- results_table[order(results_table$pval_raw),]

# The "description" attribute comes with source info, which isn't very useful, so we can get rid of this.
results_table$description <- sapply(results_table$description, function(x) str_replace(x, '[ ]*\\[Source:.*', ""))

print(head(results_table))
