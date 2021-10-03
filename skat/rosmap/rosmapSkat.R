library(SKAT)
library(RUnit)
library(VariantAnnotation)
#----------------------------------------------------------------------------------------------------
ad.snps <- c("19:45411941_T/C", "19:45412079_C/T")

khaleesi <- grepl("khaleesi", Sys.info()[["nodename"]])
hagfish <-  grepl("hagfish", Sys.info()[["nodename"]])
file.dir <- "~/github/TrenaProjectAD/explore/WGS-Harmonization/metadata-all"

specimen.file <- file.path(file.dir, "ROSMAP_biospecimen_metadata.csv")
checkTrue(file.exists(specimen.file))
tbl.spec <- read.table(specimen.file, sep=",", header=TRUE, nrow=-1)

clinical.file <- file.path(file.dir, "ROSMAP_clinical.csv")
checkTrue(file.exists(clinical.file))
tbl.clin <- read.table(clinical.file, sep=",", header=TRUE, nrow=-1)

if(khaleesi){
   dir <-  "/local/users/pshannon/github/TrenaProjectAD/explore/ampad-1898-samples/vcf"
   filename <- "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_19.recalibrated_variants.vcf.gz"
   apoe.region <- GRanges(seqnames="19", IRanges(start=45411000, end=45413000))   # 2kb
   }

if(hagfish){
   dir <-  "./"
   filename <- "../../podkat/rosmap/chr19-apoe-100k-hg19.vcf.gz"
   stopifnot(file.exists(filename))
   apoe.region <- GRanges(seqnames="19", IRanges(start=45360000, end=45460000))
   apoe.region <-
   center.base <- as.integer((45412079 + 45411941)/2)  # [1] 45412010
   shoulder <- 10000
   apoe.region <- GRanges(seqnames="19", IRanges(start=center.base-shoulder, end=center.base+shoulder))
   }

vcfFile <- file.path(dir, filename)
      # AD vcf's in hg19:
      #   rs429358   19:45411941 (GRCh37)
      #   rs7412     19:45412079 (GRCh37)

tmp.vcf <- readVcf(vcfFile, "hg19", apoe.region)
tmp.mtx.geno <- t(geno(tmp.vcf)$GT)
vcf.sampleIDs <- rownames(tmp.mtx.geno)

# 1151 genotyped samples are annotated in ROSMAP_biospecimen_metadata.csv
genotype.samples.annotated <- unique(intersect(tbl.spec$specimenID, vcf.sampleIDs))
length(genotype.samples.annotated) # [1] 1151
tbl.sampleToPatient <- subset(tbl.spec, specimenID %in% genotype.samples.annotated)[, c("specimenID", "individualID")]
checkTrue(all(tbl.sampleToPatient$specimenID %in% genotype.samples.annotated))
   # 7 patients have more than one genotype
checkEquals(length(unique(tbl.sampleToPatient$individualID)), 1144)
   # eliminate them for now
dups <- which(duplicated(tbl.sampleToPatient$individualID))
dups
tbl.sampleToPatient <- tbl.sampleToPatient[-dups, ]
checkEquals(nrow(tbl.sampleToPatient), 1144)
   # tbl.covAll will be winnowed and maybe jiggered to produce good tests below
   # by variously parameterized calls to createCovariatesTable
tbl.covAll <- merge(tbl.sampleToPatient, tbl.clin, by="individualID")
tbl.covAll$age_death[tbl.covAll$age_death == "90+"] <- "101"
tbl.covAll$age_death <- round(as.numeric(tbl.covAll$age_death), digit=0)
dim(tbl.covAll)   # 1144 19
viz <- FALSE
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_createCovariatesTable()
   test_createGeno()
   test_skatifyGenotypeMatrix()
   test_buildModelSmallerRegions()

} # runTests
#----------------------------------------------------------------------------------------------------
mapSamplePatientIDs <- function()
{

} # mapSamplePatientIDs
#----------------------------------------------------------------------------------------------------
# AD diagnosis moderately correlated to apoe4.  for testing, if injectEnrichment is TRUE,
# then align AD diagnosis to that genotype, so that podkat can find it
createCovariatesTable <- function(injectEnrichment=FALSE)
{
   coi <- c("individualID", "specimenID", "msex", "apoe_genotype", "age_death", "cogdx")
   setdiff(coi, colnames(tbl.covAll))
   #tbl.cov <- subset(tbl.covAll, Study=="ROS")
   tbl.cov <- tbl.covAll[, coi]
   dim(tbl.cov)  # 1144 6
   deleters <- which(is.na(tbl.cov$cogdx))
   length(deleters)  # 611
   if(length(deleters > 0))
        tbl.cov <- tbl.cov[-deleters,]
   dim(tbl.cov)  # 840 5

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
    checkEquals(dim(tbl.cov), c(1143, 9))
       # ad, v34 and v44 are derived columns, added by me
    coi.expected <- c("individualID", "specimenID",  "msex", "apoe_genotype", "age_death", "cogdx",
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
   vcf <- readVcf(vcfFile, "hg19", vcf.region)
   mtx.geno <- t(geno(vcf)$GT)


  #      readGenotypeMatrix(file, regions, subset,
  #                         noIndels=TRUE, onlyPass=TRUE,
  #                         na.limit=1, MAF.limit=1,
  #                         na.action=c("impute.major", "omit", "fail"),
  #                         MAF.action=c("invert", "omit", "ignore", "fail"),
  #                         sex=NULL)

   keepers <- grep("SM-", rownames(mtx.geno))
   mtx.geno <- mtx.geno[keepers,]

   name.match.indices <-  match(rownames(mtx.geno), tbl.covAll$specimenID)
   na.matches <- which(is.na(name.match.indices))
   length(na.matches)
   mtx.geno <- mtx.geno[-na.matches,]
   name.match.indices <-  match(rownames(mtx.geno), tbl.covAll$specimenID)
   stopifnot(length(which(is.na(name.match.indices))) == 0)

   rownames(mtx.geno) <- tbl.covAll$individualID[name.match.indices]
   mtx.geno

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
    checkTrue(all(grepl("19:", colnames(mtx.geno))))
    checkTrue(all(grepl("^R", rownames(mtx.geno))))
      # get just the AD-related positions, hand check the variants

      # AD vcf's in hg19:
      #   rs429358   19:45411941 (GRCh37)   ref: T
      #   rs7412     19:45412079 (GRCh37)   ref: C
    all(ad.snps %in% colnames(mtx.geno))
    mtx.adSnps <- mtx.geno[, ad.snps]
                                                          #  0/0  0/1   1/1
    checkEquals(as.numeric(table(mtx.adSnps[, ad.snps[1]])), c(847, 277,  20))
    checkEquals(as.numeric(table(mtx.adSnps[, ad.snps[2]])), c(974, 163,   7))

       # pick a double snp, make sure the patient is in tbl.cov
    pts.double <- sort(intersect(names(which(mtx.adSnps[,1] == "0/1")),
                                 names(which(mtx.adSnps[,2] == "0/1"))))
    checkEquals(length(pts.double), 21)
    pt.double.first <- head(pts.double, n=1)
    checkEquals(pt.double.first, "R2423943")
    sampleID <- subset(tbl.covAll, individualID==pt.double.first)$specimenID
    checkEquals(sampleID, "SM-CTEF9")

       # go back to the tradtional VariantAnnotation
    # checkTrue(sampleID %in% vcf.sampleIDs)
    # mtx.orig <- geno(vcf)$GT
    # khaleesi, 2kb: checkEquals(dim(mtx.orig), c(34, 1894))
    # hagfish, 100kb: checkEquals(dim(mtx.orig), c(2754, 1894))
    # roi <- unlist(lapply(ad.snps, function(loc) grep(loc, rownames(mtx.orig))))
    # checkEquals(as.character(mtx.orig[roi, sampleID]), c("0/1", "0/1"))

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
test_buildModelOneApoeRegion <- function()
{
   message(sprintf("--- test_buildModelOneApoeRegion: %d bases", width(apoe.region)))

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

} # test_buildModelOneApoeRegion
#------------------------------------------------------------------------------------------------------------------------
test_buildModelSmallerRegions <- function()
{
   roi <- GRanges(seqnames="19", IRanges(start=45408985, 45415032))
   roi <- apoe.region
   count <- 5
   count <- 20
   count <- 100
   message(sprintf("--- test_buildModelSmallerRegions: %d bases in %d blocks", width(roi), count))
   mtx.geno <- createGeno(vcf.region=roi)
   dim(mtx.geno)

   tbl.cov <- createCovariatesTable()
   shared.ids <- intersect(rownames(mtx.geno), tbl.cov$individualID)
   length(shared.ids) # 1143
   mtx.geno <- mtx.geno[shared.ids,]
   tbl.cov <- subset(tbl.cov, individualID %in% shared.ids)

   starts <- as.integer(seq(1, ncol(mtx.geno), length.out=count+1))
   ends <- starts-1
   starts <- starts[-length(starts)]
   ends <- ends[-1]
   ends[length(ends)] <- ncol(mtx.geno)
   stopifnot(sum(ends - starts + 1) == ncol(mtx.geno))

     #------------------------------------------------------------
     # gin up a strong signal, all patients with the first snp
     #------------------------------------------------------------
   has.rs429358 <- rep(0, nrow(tbl.cov))
   has.rs429358[which(mtx.geno[, "19:45411941_T/C"] != "0/0")] <- 1
   checkEquals(as.numeric(table(has.rs429358)), c(847, 296))
   tbl.cov$has.rs429358 <- has.rs429358
   null.model <- SKAT_Null_Model(has.rs429358 ~ age_death + msex, tbl.cov)
   #null.model <- SKAT_Null_Model(cogdx ~ age_death + msex, tbl.cov)

     #------------------------------------------------------------
     # another attempt at  a strong signal:
     #    only patients with cogdx==4 AND both snps
     #         patients with cogdx==1 and neither snps
     #------------------------------------------------------------
   tbl.cov.1 <-rbind(subset(tbl.cov, cogdx==4 & apoe_genotype != 33),
                     subset(tbl.cov, cogdx==1 & apoe_genotype == 33))

   null.model <- SKAT_Null_Model(cogdx ~ age_death + msex, tbl.cov.1)


   to.chrom.loc <- function(col.title){
      tokens <- strsplit(col.title, ":|_")[[1]]
      list(chrom=tokens[1], pos=as.integer(tokens[2]))
      }

   tbls.skat <- list()
   for(i in seq_len(length(starts))){
      start <- starts[i]
      end   <- ends[i]
      coi <- colnames(mtx.geno)[start:end]
      mtx.skat <- skatifyGenotypeMatrix(mtx.geno[, coi])
      mtx.skat <- mtx.skat[rownames(tbl.cov.1),]
      x <- SKAT(mtx.skat, null.model, method="SKATO")
      chrom.loc.start <- to.chrom.loc(colnames(mtx.geno)[start])
      chrom.loc.end   <- to.chrom.loc(colnames(mtx.geno)[end])
      #browser()
      score <- -log10(x$p.value)
      tbls.skat[[i]] <- data.frame(chrom=chrom.loc.start$chrom,
                                   start=chrom.loc.start$pos,
                                   end=chrom.loc.end$pos,
                                   score = score,
                                   pval = x$p.value,
                                   stringsAsFactors=FALSE)
      if(x$p.value < 0.05) printf("%d - %d: %f", start, end, x$p.value)
      }

   tbl.skat <- do.call(rbind, tbls.skat)
   tbl.skat
   dim(tbl.skat)

   if(hagfish & viz){
      igv <- start.igv("APOE", "hg19")
      tbl.adSnps <- data.frame(chrom="19",
                               start=c(45411941, 45412079),
                               end=c(45411941, 45412079),
                               stringsAsFactors=FALSE)
      track <- DataFrameAnnotationTrack("AD", tbl.adSnps, color="black")
      displayTrack(igv, track)
      #snp.counts <- colSums(mtx.geno)
      #loc.strings <- names(snp.counts)
      #tokenSet <- strsplit(loc.strings, (":|_"))
      #parseTokens <- function(tokens){
      #   data.frame(chrom=tokens[1], start=as.integer(tokens[2])-1, end=as.integer(tokens[2])+1, stringsAsFactors=FALSE)
      #   }
      #tbl.snps <- do.call(rbind, lapply(tokenSet, parseTokens))
      #tbl.snps$score <- as.integer(snp.counts)
      #track <- DataFrameQuantitativeTrack("snps", tbl.snps, color="red", autoscale=TRUE)
      #displayTrack(igv, track)

      track <- DataFrameQuantitativeTrack("SKAT",
                                          tbl.skat[, c("chrom", "start", "end", "score")], color="red", autoscale=TRUE)
      displayTrack(igv, track)
      track <- VariantTrack("rosmap", tmp.vcf)
      displayTrack(igv, track)
      #showGenomicRegion(igv, "chr19:45,411,915-45,412,103")
      }

   head(tbl.skat)
   #print(res.c)
   #plot(res.c, which="p.value")

} # test_buildModelSmallerRegions
#----------------------------------------------------------------------------------------------------
