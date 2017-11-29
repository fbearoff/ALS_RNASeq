library(DEXSeq)
library(parallel)
library(GenomicFeatures)

#define number of cores to use
BPPARAM = MulticoreParam(workers=10)

#read in sample list file
samples <- read.table(file.path(WD, "samples.txt"),
                      header = TRUE,
                      row.names = 1) 
samples$libType <- "paired-end"
samples <- samples[ which(samples$condition == "BB_90" | samples$condition == "SS17_90"), ]

#identify file locations
files <- file.path(WD, "DEXseq", samples$location)
files <- paste(files, ".txt", sep="")

#constuct DEXSeq dataset
dxd = DEXSeqDataSetFromHTSeq(countfiles = files,
                             sampleData = samples,
                             design = ~sample + exon + condition:exon,
                             flattenedfile = "~/GRCm38_construct/Mus_musculus.GRCm38.88.gff")
#normalize to read depth
dxd = estimateSizeFactors(dxd)

#estimate dispersion
dxd = estimateDispersions(dxd, BPPARAM = BPPARAM)

#infer differential exon usage
dxd = testForDEU(dxd, BPPARAM = BPPARAM)
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM = BPPARAM)

#store results
dxr = DEXSeqResults(dxd)

#generate report
mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
DEXSeqHTML(object = dxr, path="~/R_projects/ALS_RNASeq/DEXseq/BB90vsSS1790/report2",
           FDR = 0.05, color=c("red", "green"), BPPARAM = BPPARAM,
           mart = mart, filter="ensembl_gene_id",
           attributes=c("description", "mgi_symbol"))

