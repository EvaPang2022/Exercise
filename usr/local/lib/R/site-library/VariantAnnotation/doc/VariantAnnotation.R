## ----readVcF,message=FALSE----------------------------------------------------
library(VariantAnnotation)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")
vcf

## ----readVcf_showheader-------------------------------------------------------
header(vcf)

## ----headeraccessors----------------------------------------------------------
samples(header(vcf))
geno(header(vcf))

## ----readVcf_rowRanges--------------------------------------------------------
head(rowRanges(vcf), 3)

## ----readVcf_fixed------------------------------------------------------------
ref(vcf)[1:5]
qual(vcf)[1:5]

## ----readVcf_ALT--------------------------------------------------------------
alt(vcf)[1:5]

## ----geno_hdr-----------------------------------------------------------------
geno(vcf)
sapply(geno(vcf), class)

## ----explore_geno-------------------------------------------------------------
geno(header(vcf))["DS",]

## ----dim_geno-----------------------------------------------------------------
DS <-geno(vcf)$DS
dim(DS)
DS[1:3,]

## ----fivenum------------------------------------------------------------------
fivenum(DS)

## ----DS_zero------------------------------------------------------------------
length(which(DS==0))/length(DS)

## ----DS_hist, fig=TRUE--------------------------------------------------------
hist(DS[DS != 0], breaks=seq(0, 2, by=0.05),
    main="DS non-zero values", xlab="DS")

## ----info---------------------------------------------------------------------
info(vcf)[1:4, 1:5]

## ----examine_dbSNP, message=FALSE, warning=FALSE------------------------------
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
vcf_rsids <- names(rowRanges(vcf))
chr22snps <- snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, "22")
chr22_rsids <- mcols(chr22snps)$RefSNP_id
in_dbSNP <- vcf_rsids %in% chr22_rsids
table(in_dbSNP) 

## ----header_info--------------------------------------------------------------
info(header(vcf))[c("VT", "LDAF", "RSQ"),]

## ----examine_quality----------------------------------------------------------
metrics <- data.frame(QUAL=qual(vcf), in_dbSNP=in_dbSNP,
    VT=info(vcf)$VT, LDAF=info(vcf)$LDAF, RSQ=info(vcf)$RSQ)

## ----examine_ggplot2, message=FALSE, warning=FALSE, fig=TRUE------------------
library(ggplot2)
ggplot(metrics, aes(x=RSQ, fill=in_dbSNP)) +
    geom_density(alpha=0.5) +
    scale_x_continuous(name="MaCH / Thunder Imputation Quality") +
    scale_y_continuous(name="Density") +
    theme(legend.position="top")

## ----subset_ranges------------------------------------------------------------
rng <- GRanges(seqnames="22", ranges=IRanges(
           start=c(50301422, 50989541), 
           end=c(50312106, 51001328),
           names=c("gene_79087", "gene_644186")))

## ----subset_TabixFile---------------------------------------------------------
tab <- TabixFile(fl)
vcf_rng <- readVcf(tab, "hg19", param=rng)

## -----------------------------------------------------------------------------
head(rowRanges(vcf_rng), 3)

## ----subset_scanVcfHeader-----------------------------------------------------
hdr <- scanVcfHeader(fl)
## e.g., INFO and GENO fields
head(info(hdr), 3)
head(geno(hdr), 3)

## ----subset_ScanVcfParam------------------------------------------------------
## Return all 'fixed' fields, "LAF" from 'info' and "GT" from 'geno'
svp <- ScanVcfParam(info="LDAF", geno="GT")
vcf1 <- readVcf(fl, "hg19", svp)
names(geno(vcf1))

## ----subset_ScanVcfParam_new--------------------------------------------------
svp_all <- ScanVcfParam(info="LDAF", geno="GT", which=rng)
svp_all

## ----locate_rename_seqlevels, message=FALSE, warning=FALSE--------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(vcf) <- "chr22"
rd <- rowRanges(vcf)
loc <- locateVariants(rd, txdb, CodingVariants())
head(loc, 3)

## ----AllVariants, eval=FALSE--------------------------------------------------
#  allvar <- locateVariants(rd, txdb, AllVariants())

## ----locate_gene_centric------------------------------------------------------
## Did any coding variants match more than one gene?
splt <- split(mcols(loc)$GENEID, mcols(loc)$QUERYID) 
table(sapply(splt, function(x) length(unique(x)) > 1))

## Summarize the number of coding variants by gene ID.
splt <- split(mcols(loc)$QUERYID, mcols(loc)$GENEID)
head(sapply(splt, function(x) length(unique(x))), 3)

## ----predictCoding, warning=FALSE---------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)
coding[5:7]

## ----predictCoding_frameshift-------------------------------------------------
## CONSEQUENCE is 'frameshift' where translation is not possible
coding[mcols(coding)$CONSEQUENCE == "frameshift"]

## ----nonsynonymous------------------------------------------------------------
nms <- names(coding) 
idx <- mcols(coding)$CONSEQUENCE == "nonsynonymous"
nonsyn <- coding[idx]
names(nonsyn) <- nms[idx]
rsids <- unique(names(nonsyn)[grep("rs", names(nonsyn), fixed=TRUE)])

## ----polyphen, message=FALSE, warning=FALSE-----------------------------------
library(PolyPhen.Hsapiens.dbSNP131)
pp <- select(PolyPhen.Hsapiens.dbSNP131, keys=rsids,
          cols=c("TRAININGSET", "PREDICTION", "PPH2PROB"))
head(pp[!is.na(pp$PREDICTION), ]) 

## ----snpMatrix, message=FALSE-------------------------------------------------
res <- genotypeToSnpMatrix(vcf) 
res

## ----snpMatrix_ALT------------------------------------------------------------
allele2 <- res$map[["allele.2"]]
## number of alternate alleles per variant
unique(elementNROWS(allele2))

## ----message=FALSE------------------------------------------------------------
fl.gl <- system.file("extdata", "gl_chr1.vcf", package="VariantAnnotation")
vcf.gl <- readVcf(fl.gl, "hg19")
geno(vcf.gl)

## Convert the "GL" FORMAT field to a SnpMatrix
res <- genotypeToSnpMatrix(vcf.gl, uncertain=TRUE)
res
t(as(res$genotype, "character"))[c(1,3,7), 1:5]

## Compare to a SnpMatrix created from the "GT" field
res.gt <- genotypeToSnpMatrix(vcf.gl, uncertain=FALSE)
t(as(res.gt$genotype, "character"))[c(1,3,7), 1:5]

## What are the original likelihoods for rs58108140?
geno(vcf.gl)$GL["rs58108140", 1:5]

## ----writeVcf, message=FALSE, warning=FALSE-----------------------------------
fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
out1.vcf <- tempfile()
out2.vcf <- tempfile()
in1 <- readVcf(fl, "hg19")
writeVcf(in1, out1.vcf)
in2 <- readVcf(out1.vcf, "hg19")
writeVcf(in2, out2.vcf)
in3 <- readVcf(out2.vcf, "hg19")
identical(rowRanges(in1), rowRanges(in3))
identical(geno(in1), geno(in2))

## ----eval=FALSE---------------------------------------------------------------
#  readVcf(TabixFile(fl, yieldSize=10000))

## ----eval=FALSE---------------------------------------------------------------
#  readVcf(TabixFile(fl), param=ScanVcfParam(info='DP', geno='GT'))

## ----eval=FALSE---------------------------------------------------------------
#  readGT(fl)

## ----eval=FALSE---------------------------------------------------------------
#  library(microbenchmark)
#  fl <- "ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
#  ys <- c(100, 1000, 10000, 100000)
#  
#  ## readGT() input only 'GT':
#  fun <- function(fl, yieldSize) readGT(TabixFile(fl, yieldSize))
#  lapply(ys, function(i) microbenchmark(fun(fl, i), times=5)
#  
#  ## readVcf() input only 'GT' and 'ALT':
#  fun <- function(fl, yieldSize, param)
#             readVcf(TabixFile(fl, yieldSize), "hg19", param=param)
#  param <- ScanVcfParam(info=NA, geno="GT", fixed="ALT")
#  lapply(ys, function(i) microbenchmark(fun(fl, i, param), times=5)
#  
#  ## readVcf() input all variables:
#  fun <- function(fl, yieldSize) readVcf(TabixFile(fl, yieldSize), "hg19")
#  lapply(ys, function(i) microbenchmark(fun(fl, i), times=5))

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

