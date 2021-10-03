library(podkat)
library(RUnit)
library(VariantAnnotation)
#----------------------------------------------------------------------------------------------------
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
   filename <- "chr19-apoe-100k-hg19.vcf.gz"
   apoe.region <- GRanges(seqnames="19", IRanges(start=45360000, end=45460000))
   apoe.region <-
   center.base <- as.integer((45412079 + 45411941)/2)  # [1] 45412010
   apoe.region <- GRanges(seqnames="19", IRanges(start=center.base-1000, end=center.base+1000))
   }

vcfFile <- file.path(dir, filename)
      # AD vcf's in hg19:
      #   rs429358   19:45411941 (GRCh37)
      #   rs7412     19:45412079 (GRCh37)

vcf <- readVcf(vcfFile, "hg19", apoe.region)
vcf.sampleIDs <- colnames(geno(vcf)$GT)

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
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_createCovariatesTable()
   test_createGeno()

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

   has.rs429358 <- rep(0, nrow(geno))
   has.rs429358[(which(as.matrix(geno)[, "19:45411941_T/C"] > 0))] <- 1

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
createGeno <- function()
{
  #      readGenotypeMatrix(file, regions, subset,
  #                         noIndels=TRUE, onlyPass=TRUE,
  #                         na.limit=1, MAF.limit=1,
  #                         na.action=c("impute.major", "omit", "fail"),
  #                         MAF.action=c("invert", "omit", "ignore", "fail"),
  #                         sex=NULL)

   geno <- readGenotypeMatrix(vcfFile, apoe.region)
   keepers <- grep("SM-", rownames(geno))
   geno <- geno[keepers,]

   name.match.indices <-  match(rownames(geno), tbl.covAll$specimenID)
   na.matches <- which(is.na(name.match.indices))
   length(na.matches)
   geno <- geno[-na.matches,]
   name.match.indices <-  match(rownames(geno), tbl.covAll$specimenID)
   stopifnot(length(which(is.na(name.match.indices))) == 0)

   rownames(geno) <- tbl.covAll$individualID[name.match.indices]
   geno

} # createGeno
#----------------------------------------------------------------------------------------------------
test_createGeno <- function()
{
    message(sprintf("--- test_createGeno"))

    geno <- createGeno()
    # khaleesi, 2kb checkEquals(dim(geno), c(1144, 32))
    # hagfish, 100kb, 2288
    colnames(geno)
    checkTrue(all(grepl("19:", colnames(geno))))
    checkTrue(all(grepl("^R", rownames(geno))))
      # get just the AD-related positions, hand check the variants

      # AD vcf's in hg19:
      #   rs429358   19:45411941 (GRCh37)   ref: T
      #   rs7412     19:45412079 (GRCh37)   ref: C
    ad.snps <- c("19:45411941_T/C", "19:45412079_C/T")
    all(ad.snps %in% colnames(geno))
    mtx.geno <- as.matrix(geno)[, ad.snps]
    sums.by.row <- rowSums(mtx.geno)
    table(rowSums(mtx.geno))
       #   0   1   2
       # 698 398  48
    checkEquals(length(which(sums.by.row==2)), 48)
       # choose the first sample with variants in both positions
    x <- names(which(sums.by.row == 2)[1])
    checkEquals(x, "R6253512")
    sampleID <- subset(tbl.covAll, individualID==x)$specimenID
    checkEquals(sampleID, "SM-CJEGR")

       # go back to the tradtional VariantAnnotation
    checkTrue(sampleID %in% vcf.sampleIDs)
    mtx.orig <- geno(vcf)$GT
    # khaleesi, 2kb: checkEquals(dim(mtx.orig), c(34, 1894))
    # hagfish, 100kb: checkEquals(dim(mtx.orig), c(2754, 1894))
    roi <- unlist(lapply(ad.snps, function(loc) grep(loc, rownames(mtx.orig))))
    checkEquals(as.character(mtx.orig[roi, sampleID]), c("0/1", "0/1"))

} # test_createGeno
#----------------------------------------------------------------------------------------------------
test_buildModel <- function()
{
   tbl.cov <- createCovariatesTable()
   geno <- createGeno()
   dim(geno)

   shared.ids <- intersect(rownames(geno), tbl.cov$individualID)
   length(shared.ids) # 1143
   geno <- geno[shared.ids,]
   tbl.cov <- subset(tbl.cov, individualID %in% shared.ids)
   checkEquals(nrow(tbl.cov), 1143)
   checkEquals(nrow(geno), 1143)

   null.model <- nullModel(has.rs429358 ~ age_death + msex, tbl.cov)

   lm.model <- lm(has.rs429358 ~ age_death + msex, tbl.cov)
   summary(lm.model)  # adjusted R-squared: useless model, as expected  0.008743

   #gr.oi <- GRanges(seqnames="19", IRanges(start=45407983, end=45416174))
   gr.oi <- apoe.region
   windows <- partitionRegions(gr.oi, width=250, overlap=0.5)

   res.c <- assocTest(geno, null.model, windows)
   tbl.res <- as.data.frame(res.c)
   tbl.res$score <- -log10(tbl.res$p.value)
   tbl.res$seqnames <- as.character(tbl.res$seqnames)
   fivenum(tbl.res$score)
   tbl.res

   if(hagfish){
      igv <- start.igv("APOE", "hg19")
      tbl.adSnps <- data.frame(chrom="19",
                               start=c(45411941, 45412079),
                               end=c(45411941, 45412079),
                               stringsAsFactors=FALSE)
      track <- DataFrameAnnotationTrack("AD", tbl.adSnps, color="black")
      displayTrack(igv, track)
      snp.counts <- colSums(geno)
      loc.strings <- names(snp.counts)
      tokenSet <- strsplit(loc.strings, (":|_"))
      parseTokens <- function(tokens){
         data.frame(chrom=tokens[1], start=as.integer(tokens[2])-1, end=as.integer(tokens[2])+1, stringsAsFactors=FALSE)
         }
      tbl.snps <- do.call(rbind, lapply(tokenSet, parseTokens))
      tbl.snps$score <- as.integer(snp.counts)
      track <- DataFrameQuantitativeTrack("snps", tbl.snps, color="red", autoscale=TRUE)
      displayTrack(igv, track)

      track <- DataFrameQuantitativeTrack("podkat", tbl.res[, c("seqnames", "start", "end", "score")], color="red", autoscale=TRUE)
      displayTrack(igv, track)
      track <- VariantTrack("rosmap", vcf)
      displayTrack(igv, track)
      showGenomicRegion(igv, "chr19:45,411,915-45,412,103")
      }

   head(tbl.res)
   print(res.c)
   plot(res.c, which="p.value")

} # buildModel
#----------------------------------------------------------------------------------------------------
