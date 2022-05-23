setwd("/scratch/ahc87874/Check/genoQC")
library(gds2bgen)

for (i in 1:22) {

bgen_fn <- paste("chr", i, ".bgen", sep = "")
gds_out <- paste("chr", i, ".gds", sep = "")    
#seqBGEN_Info(bgen_fn)

seqBGEN2GDS(bgen_fn, gds_out,
    storage.option="ZIP_RA",  # compression option, e.g., ZIP_RA for zlib or LZ4_RA for LZ4
    float.type="packed8",      # 8-bit packed real numbers
    geno=TRUE,     # 2-bit integer genotypes, stored in 'genotype/data'
    dosage=TRUE,    # numeric alternative allele dosages, stored in 'annotation/format/DS'
    prob=FALSE,     # numeric genotype probabilities, stored in 'annotation/format/GP'
    parallel=4      # the number of cores
)

#library(SeqArray)
#(f <- seqOpen("chr22.gds"))
#seqClose(f)
}
