library(EndophenotypeExplorer)
library(SKAT)
library(RUnit)
if(!exists("etx"))
    etx <- EndophenotypeExplorer$new("NDUFS2", "hg38", initialize.snpLocs=TRUE)

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_createCovariatesTable()
   test_createGeno()

} # runTests
#----------------------------------------------------------------------------------------------------
# AD diagnosis moderately correlated to apoe4.  for testing, if injectEnrichment is TRUE,
# then align AD diagnosis to that genotype, so that podkat can find it
createCovariatesTable <- function(injectEnrichment=FALSE)
{
   tbl.pt <- etx$get.rosmap.patient.table(NA)
   coi <- c("individualID", "msex", "apoe_genotype", "age_death", "cogdx")
   #setdiff(coi, colnames(tbl.covAll))
   #tbl.cov <- subset(tbl.covAll, Study=="ROS")
   tbl.cov <- tbl.pt[, coi]
   dim(tbl.cov)  # 3583 5
   deleters <- which(is.na(tbl.cov$cogdx))
   length(deleters)  # 1762
   if(length(deleters > 0))
        tbl.cov <- tbl.cov[-deleters,]
   dim(tbl.cov)  # 1821 5

     # choose a first, simple outcomes variable. cogdx?
     #    1     NCI, No cognitive impairment (No impaired domains)
     #    2     MCI, Mild cognitive impairment (One impaired domain) and NO other cause of CI
     #    3     MCI, Mild cognitive impairment (One impaired domain) AND another cause of CI
     #    4     AD, Alzheimer's disease and NO other cause of CI (NINCDS PROB AD)
     #    5     AD, Alzheimer's disease AND another cause of CI (NINCDS POSS AD)
     #    6     Other dementia. Other primary cause of dementia
     #
     # reduce these to 0 1 1 5 4 2xo

   outcomes <- tbl.cov$cogdx
   tbl.cov$ad <- rep(0, nrow(tbl.cov))
   ad.patients <- c(which(tbl.cov$cogdx==4), which(tbl.cov$cogdx==5))
   length(ad.patients)
   nrow(tbl.cov)
   tbl.cov$ad[ad.patients] <- 1
   tbl.cov$v34 <- rep(0, nrow(tbl.cov))
   tbl.cov$v34[which(tbl.cov$apoe_genotype==34)] <- 1
   tbl.cov$v44 <- rep(0, nrow(tbl.cov))
   tbl.cov$v44[which(tbl.cov$apoe_genotype==44)] <- 1

     # fabricate an outcomes variable from actual variants at rs29358
     #   rs429358   19:45411941 (GRCh37)
     #   rs7412     19:45412079 (GRCh37)


   rownames(tbl.cov) <- tbl.cov$individualID

   tbl.cov

} # createCovariatesTable
#----------------------------------------------------------------------------------------------------
test_createCovariatesTable <- function()
{
    message(sprintf("--- test_createCovariatesTable"))
    tbl.cov <- createCovariatesTable(injectEnrichment=TRUE)
    checkEquals(dim(tbl.cov), c(1821, 8))
       # ad, v34 and v44 are derived columns, added by me
    coi.expected <- c("individualID", "msex", "apoe_genotype", "age_death", "cogdx",
                      "ad", "v34", "v44")
    checkEquals(colnames(tbl.cov), coi.expected)
    checkTrue(all(grepl("^R", rownames(tbl.cov))))

       # v44 has 11x increased AD risk,  do we see that?

    checkEquals(median(subset(tbl.cov, v44==1)$cogdx), 4)  # AD, no other causes of cog impairment
    checkEquals(median(subset(tbl.cov, v44==0)$cogdx), 2)  # mild cognitive impairment

} # test_createCovariatesTable
#----------------------------------------------------------------------------------------------------
createGeno <- function(vcf.region=GRanges(seqnames="19", IRanges(start=45411000, end=45413000)))
{
      # get mtx.geno for all rsids with pval < 0.0001
   pval.cutoff <- 1e-6
   tbl.eqtl <- etx$getEQTLsForGene()
   tbl.eqtl <- subset(tbl.eqtl, pvalue <= pval.cutoff & study=="ampad-rosmap")
   dim(tbl.eqtl)  # 19 10
   rsids <- tbl.eqtl$rsid
   mtx.geno <- etx$getGenoMatrixByRSID(rsids)
   dim(mtx.geno)   # 19 1894
   rsid.rownames <- etx$locsToRSID(rownames(mtx.geno), "hg19")
   rownames(mtx.geno) <- rsid.rownames
   tbl.map <- etx$getIdMap()
     # convert vcf sample names to corresponding rosmap patient names
   tbl.rosmap.vcf <- subset(tbl.map, study=="rosmap" & assay=="vcf")
   rosmap.vcf.samples <- intersect(tbl.rosmap.vcf$sample, colnames(mtx.geno))  # 1151
   mtx.geno <- mtx.geno[, rosmap.vcf.samples]
   indices <- match(tbl.rosmap.vcf$sample, colnames(mtx.geno))
   length(indices)
   colnames(mtx.geno) <- tbl.rosmap.vcf$patient[indices]
   mtx.geno.t <- t(mtx.geno)
   invisible(mtx.geno.t)

} # createGeno
#----------------------------------------------------------------------------------------------------
test_createGeno <- function()
{
    message(sprintf("--- test_createGeno"))

    mtx.geno <- createGeno()
    checkTrue(is.matrix(mtx.geno))
    # khaleesi, 2kb checkEquals(dim(mtx.geno), c(1144, 32))
    # hagfish, 100kb, 2288
    colnames(mtx.geno)
    checkTrue(all(grepl("^rs", colnames(mtx.geno))))
    checkTrue(all(grepl("^R", rownames(mtx.geno))))
      # get just the AD-related positions, hand check the variants

      # AD vcf's in hg19:
      #   rs429358   19:45411941 (GRCh37)   ref: T
      #   rs7412     19:45412079 (GRCh37)   ref: C

    checkEquals(as.numeric(table(mtx.geno[, 1])), c(1016, 131, 4))

       # hand-check some of the eQTLS
    rsid.strongest.eqtl <- "rs1136224"   # has pvalue of 1e-18
    rsid.strongest.gwas <- "rs4575098"

        # counts of 0/0, 0/1, 1/1
    checkEquals(as.numeric(table(mtx.geno[, rsid.strongest.eqtl])), c(827, 296, 28))
    checkEquals(as.numeric(table(mtx.geno[, rsid.strongest.gwas])), c(686, 405, 60))

    mtx.geno.test <- etx$getGenoMatrixByRSID(c(rsid.strongest.eqtl, rsid.strongest.gwas))
    new.names <- etx$locsToRSID(rownames(mtx.geno.test), "hg19")
        # mtx.geno.test is row-sorted by ascending chrom loc
    checkEquals(as.character(new.names), c(rsid.strongest.gwas, rsid.strongest.eqtl))
    rownames(mtx.geno.test) <- as.character(new.names)

    vcf.samples <- colnames(mtx.geno.test)
    rosmap.vcf.samples <- intersect(vcf.samples, subset(tbl.map, study=="rosmap" & assay=="vcf")$sample)
    checkEquals(length(rosmap.vcf.samples), 1151)
    mtx.geno.test <- mtx.geno.test[, rosmap.vcf.samples]

        # a fresh extraction, and mtx.geno: both should have 60 homozygous patients
    rs4575098.homs.test <- which(mtx.geno.test["rs4575098",] == "1/1")
    checkEquals(length(rs4575098.homs.test), 60)
    checkEquals(length(which(mtx.geno[, "rs4575098"] == "1/1")), 60)

        # just 28 for the most significant eQTL.  are all "1/1" in mtx.geno and mtx.geno.test?
    rs1136224.homs.test <- which(mtx.geno.test["rs1136224",] == "1/1")
    checkEquals(length(rs1136224.homs.test), 28)
    homozygous.patients <- names(which(mtx.geno[, "rs1136224"] == "1/1"))
    checkEquals(length(homozygous.patients), 28)
    homozygous.vcf.sample.names <- subset(tbl.map, patient %in% homozygous.patients &
                                                   assay=="vcf" & study=="rosmap")$sample
    checkEquals(length(homozygous.vcf.sample.names), 28)
    checkTrue(all(mtx.geno.test["rs1136224", homozygous.vcf.sample.names] == "1/1"))

} # test_createGeno
#----------------------------------------------------------------------------------------------------
# create a numeric genotype matrix with each row as a different individual and each column
# as a separate gene/snp.  Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA,
# Aa, aa, and missing, where A is a major allele and a is a minor allele.  Missing
# genotypes will be imputed by the simple Hardy-Weinberg equilibrium (HWE) based
# imputation.
skatifyGenotypeMatrix <- function(mtx.raw)
{
   #stopifnot(all(rownames(mtx.raw) %in% tbl.cov$individualID))

   noReads <- which(mtx.raw == "./.")
   AA <- which(mtx.raw == "0/0")
   aa <- which(mtx.raw == "1/1")
   Aa <- grep("0/[1-9]", mtx.raw)
   head(setdiff(seq_len(length(mtx.raw)), c(AA, noReads, aa, Aa)))
      # see this discussion of unusual genotypes in rosmap vcf file
      # ~/github/notes/log
      # --- skat vs actual genotype from rosmap: 98% fit, 2% do not, as
      #  discovered in ~/github/SKAT/inst/unitTests/test_rosmapSkat.R

   #stopifnot(sum(length(AA),
   #              length(noReads),
   #              length(aa),
   #              length(Aa)) > 0.98 * length(mtx.raw))
   # wt <- c(hets, noReads)
   # variants <- setdiff(seq_len(length(mtx.raw)), wt)

   mtx.out <- matrix(0, nrow=nrow(mtx.raw), ncol=ncol(mtx.raw))
   mtx.out[noReads] <- 9
   mtx.out[AA] <- 0
   mtx.out[Aa] <- 1
   mtx.out[aa] <- 2
   rownames(mtx.out) <- rownames(mtx.raw)
   colnames(mtx.out) <- colnames(mtx.raw)

   #patientIDs <- tbl.cov$patientID[match(colnames(mtx.out), tbl.cov$specimen)]
   #browser()
   #stopifnot(length(patientIDs) == ncol(mtx.out))
   #colnames(mtx.out) <- patientIDs

   invisible(mtx.out)

} # skatifyGenotypeMatrix
#------------------------------------------------------------------------------------------------------------------------
test_skatifyGenotypeMatrix <- function()
{
   message(sprintf("--- test_skatifyGenotypeMatrix"))

   mtx.geno <- createGeno()
   dim(mtx.geno)
   mtx.test <- mtx.geno[2:5, c("19:45411934_G/A", "19:45411941_T/C", "19:45411987_G/A", "19:45412040_C/T")]
   mtx.test [2,3] <- "./."   # inject a no-reads value
   mtx.test [3,2] <- "1/1"   # inject a homozygous variant value

                                             # ./. 0/0 0/1 1/1
   checkEquals(as.numeric(table(mtx.test)), c(1, 13,  1, 1))

   mtx <- skatifyGenotypeMatrix(mtx.test)
                                          #  0  1  2  9
   checkEquals(as.numeric(table(mtx)), c(13, 1, 1, 1))

   checkEquals(dim(mtx), dim(mtx.test))
   checkEquals(rownames(mtx), c("R6478102", "R4567280", "R1980623", "R8444624"))
   checkEquals(colnames(mtx), c("19:45411934_G/A", "19:45411941_T/C", "19:45411987_G/A", "19:45412040_C/T"))

   mtx.test <- mtx.geno[1000:1020, 1:34]
   checkEquals(dim(mtx.test), c(21, 34))
   mtx <- skatifyGenotypeMatrix(mtx.test)
   checkEquals(dim(mtx), dim(mtx.test))
   checkEquals(as.numeric(table(mtx)), c(693, 21)) # AA, Aa

      # now the entire mtx.geno
   mtx <- skatifyGenotypeMatrix(mtx.geno)
   checkEquals(dim(mtx), dim(mtx.geno))
                                      #     0    1   2
   checkEquals(as.numeric(table(mtx)), c(38109, 756, 31))

   checkEquals(head(colnames(mtx), n=3), c("19:45411042_G/A", "19:45411045_G/A", "19:45411110_T/C"))
   checkEquals(head(rownames(mtx), n=3), c("R1977848", "R6478102", "R4567280"))

} # test_skatifyGenotypeMatrix
#------------------------------------------------------------------------------------------------------------------------
ndufs2 <- function()
{
   tbl.cov <- createCovariatesTable()
   mtx.geno <- createGeno()
   dim(mtx.geno)

   shared.ids <- intersect(rownames(mtx.geno), tbl.cov$individualID)
   length(shared.ids) # 1143
   mtx.geno <- mtx.geno[shared.ids,]
   tbl.cov <- subset(tbl.cov, individualID %in% shared.ids)
   checkEquals(nrow(tbl.cov), 1143)
   checkEquals(nrow(mtx.geno), 1143)

     # gin up a strong signal, all patients with the first snp
   has.rs429358 <- rep(0, nrow(tbl.cov))
   has.rs429358[which(mtx.geno[, "19:45411941_T/C"] != "0/0")] <- 1
   checkEquals(as.numeric(table(has.rs429358)), c(847, 296))
   tbl.cov$has.rs429358 <- has.rs429358

   null.model <- SKAT_Null_Model(has.rs429358 ~ age_death + msex, tbl.cov)
   lm.model <- lm(has.rs429358 ~ age_death + msex, tbl.cov)
   checkTrue(abs(summary(lm.model)$adj.r.squared) < 0.01)   # no variance accounted for
   summary(lm.model)  # adjusted R-squared: useless model, as expected  0.008743

   mtx.skat <- skatifyGenotypeMatrix(mtx.geno)
   x <- SKAT(mtx.skat, null.model)
   checkTrue(x$p.value < 1e-60)
   x.0 <- SKAT(mtx.skat, null.model, method="SKATO")
   x.0$p.value

} # dufs2
#----------------------------------------------------------------------------------------------------

