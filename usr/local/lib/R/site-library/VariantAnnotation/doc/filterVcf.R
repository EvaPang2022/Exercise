## ----eval=FALSE---------------------------------------------------------------
#  tabix.file <- TabixFile(file.gz, yieldSize=10000)
#  filterVcf(tabix.file, genome, destination.file,
#                prefilters=prefilters, filters=filters)

## ----prefilters---------------------------------------------------------------
isGermlinePrefilter <- function(x) {
    grepl("Germline", x, fixed=TRUE)
}

notInDbsnpPrefilter <- function(x) {
    !(grepl("dbsnp", x, fixed=TRUE))
}

## ----filters------------------------------------------------------------------
## We will use isSNV() to filter only SNVs

allelicDepth <- function(x) {
    ##  ratio of AD of the "alternate allele" for the tumor sample
    ##  OR "reference allele" for normal samples to total reads for
    ##  the sample should be greater than some threshold (say 0.1,
    ##  that is: at least 10% of the sample should have the allele
    ##  of interest)
    ad <- geno(x)$AD
    tumorPct <- ad[,1,2,drop=FALSE] / rowSums(ad[,1,,drop=FALSE])
    normPct <- ad[,2,1, drop=FALSE] / rowSums(ad[,2,,drop=FALSE])
    test <- (tumorPct > 0.1) | (normPct > 0.1)
    as.vector(!is.na(test) & test)
}

## ----createFilterRules, message=FALSE, warning=FALSE--------------------------
library(VariantAnnotation)
prefilters <- FilterRules(list(germline=isGermlinePrefilter,
                               dbsnp=notInDbsnpPrefilter))
filters <- FilterRules(list(isSNV=isSNV, AD=allelicDepth))

## ----createFilteredFile, message=FALSE, warning=FALSE-------------------------
file.gz     <- system.file("extdata", "chr7-sub.vcf.gz",
                           package="VariantAnnotation")
file.gz.tbi <- system.file("extdata", "chr7-sub.vcf.gz.tbi",
                           package="VariantAnnotation")
destination.file <- tempfile()
tabix.file <- TabixFile(file.gz, yieldSize=10000)
filterVcf(tabix.file,  "hg19", destination.file,
          prefilters=prefilters, filters=filters, verbose=TRUE)

## ----mcf7regulatoryRegions, message=FALSE, warning=FALSE----------------------
library(AnnotationHub)
hub <- AnnotationHub()
id <- names(query(hub, "wgEncodeUwTfbsMcf7CtcfStdPkRep1.narrowPeak"))
mcf7.gr <- hub[[tail(id, 1)]]

## ----findOverlaps-------------------------------------------------------------
vcf <- readVcf(destination.file, "hg19")
seqlevels(vcf) <- paste("chr", seqlevels(vcf), sep="")
ov.mcf7 <- findOverlaps(vcf, mcf7.gr)

## ----locateVariant, message=FALSE, warning=FALSE------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
locateVariants(vcf[6,], txdb, AllVariants())

## ----eval=FALSE---------------------------------------------------------------
#  library(VariantAnnotation)
#  file.gz <- "somaticVcfBeta-HCC1187-H-200-37-ASM-T1-N1.vcf.gz"
#  stopifnot(file.exists(file.gz))
#  file.gz.tbi <- paste(file.gz, ".tbi", sep="")
#  if(!(file.exists(file.gz.tbi)))
#      indexTabix(file.gz, format="vcf")
#  start.loc <- 55000000
#  end.loc   <- 56000000
#  chr7.gr <- GRanges("7", IRanges(start.loc, end.loc))
#  params <- ScanVcfParam(which=chr7.gr)
#  vcf <- readVcf(TabixFile(file.gz), "hg19", params)
#  writeVcf(vcf, "chr7-sub.vcf")
#  bgzip("chr7-sub.vcf", overwrite=TRUE)
#  indexTabix("chr7-sub.vcf.gz", format="vcf")

