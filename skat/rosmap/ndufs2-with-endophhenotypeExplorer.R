library(EndophenotypeExplorer)
library(SKAT)
library(RUnit)
library(plyr)
library(logicFS)
if(!exists("my.anneal"))
   my.anneal <- logreg.anneal.control(start = 2, end = -2, iter = 10000)

if(!exists("etx")){
    etx <- EndophenotypeExplorer$new("NDUFS2", "hg38", initialize.snpLocs=TRUE)
    tbl.eqtl <- etx$getEQTLsForGene()
    }
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_createCovariatesTable()
   test_createGeno()

} # runTests
#----------------------------------------------------------------------------------------------------
#  For the ROSMAP study, we defined
#  AD cases:  individuals with a Braak greater than or equal to 4,
#             CERAD score less than or equal to 2,
#              and a cognitive diagnosis of probable AD with no other causes (cogdx = 4)
#  controls: individuals with Braak less than or  equal to 3, CERAD score greater
#  than or equal to 3, and cognitive diagnosis of ‘no cognitive impairment’ (cogdx = 1).


createCovariatesTable <- function(injectEnrichment=FALSE)
{
   tbl.pt <- etx$get.rosmap.patient.table(NA)
   tbl.pt$ad  <- with(tbl.pt, braaksc >= 4 & ceradsc <= 2 & cogdx == 4)
   tbl.pt$ctl <- with(tbl.pt, braaksc <= 3 & ceradsc >= 3 & cogdx == 1)
   tbl.pt$dx <- -1
   tbl.pt$dx[tbl.pt$ad] <- 1
   tbl.pt$dx[tbl.pt$ctl] <- 0
   table(tbl.pt$dx)
     #   -1  0  1
     #  2935   203 445

   coi <- c("individualID", "msex", "apoe_genotype", "age_death", "cogdx", "ad", "ctl", "dx",
            "ceradsc", "braaksc")
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

   no.clear.diagnosis <- which(tbl.cov$dx < 0)
   dim(tbl.cov)
   length(no.clear.diagnosis)
   tbl.cov <- tbl.cov[-no.clear.diagnosis,]
   tbl.cov

} # createCovariatesTable
#----------------------------------------------------------------------------------------------------
test_createCovariatesTable <- function()
{
    message(sprintf("--- test_createCovariatesTable"))
    tbl.cov <- createCovariatesTable(injectEnrichment=TRUE)
    checkEquals(dim(tbl.cov), c(648, 12))
    #checkEquals(dim(tbl.cov), c(1821, 10))
       # ad, v34 and v44 are derived columns,  added by me
    coi.expected <- c("individualID","msex","apoe_genotype","age_death","cogdx","ad","ctl","dx","v34","v44", "dx")
    checkTrue(all(coi.expected %in% colnames(tbl.cov)))
    checkTrue(all(grepl("^R", rownames(tbl.cov))))

       # v44 has 11x increased AD risk,  do we see that?

    #checkEquals(median(subset(tbl.cov, v44==1)$cogdx), 4)  # AD, no other causes of cog impairment
    #checkEquals(median(subset(tbl.cov, v44==0)$cogdx), 2)  # mild cognitive impairment

} # test_createCovariatesTable
#----------------------------------------------------------------------------------------------------
createGeno <- function(pval.cutoff=1e-6)
{
   tbl.eqtl <- subset(tbl.eqtl, pvalue <= pval.cutoff & study=="ampad-rosmap")
   dim(tbl.eqtl)  # 19 10
   rsids <- tbl.eqtl$rsid
   mtx.geno <- etx$getGenoMatrixByRSID(rsids)
   dim(mtx.geno)   # 19 1894
   rsid.rownames <- etx$locsToRSID(rownames(mtx.geno), "hg19")
   failures <- which(is.na(names(rsid.rownames)))
   if(length(failures) > 0)
       rsid.rownames <- rsid.rownames[-failures]
   if(length(rsid.rownames) != nrow(mtx.geno))
       browser()
   rownames(mtx.geno) <- as.character(rsid.rownames)
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
    tbl.map <- etx$getIdMap()
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

    save(mtx.geno, file=sprintf("mtx.geno-%s.RData", gsub(" ", ".", Sys.time(), fixed=TRUE)))

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

   mtx.geno <- get(load("mtx.geno-2021-10-03.17:18:06.RData"))
   dim(mtx.geno)   # 1151 19
   mtx.test <- mtx.geno[2:5, 1:8]
   mtx.test [2,3] <- "./."   # inject a no-reads value
   mtx.test [3,2] <- "1/1"   # inject a homozygous variant value

                                             # ./. 0/0 0/1 1/1
   checkEquals(as.numeric(table(mtx.test)), c(1, 26, 3, 2))

   mtx <- skatifyGenotypeMatrix(mtx.test)
                                       #  0  1  2  9
   checkEquals(as.numeric(table(mtx)), c(26, 3, 2, 1))

   checkEquals(dim(mtx), dim(mtx.test))
   checkEquals(rownames(mtx), c("R6478102", "R4567280", "R1980623", "R8444624"))
   checkEquals(colnames(mtx), c("rs17356051", "rs145282062", "rs4575098", "rs11585858",
                                "rs4233366", "rs34884997", "rs35630547", "rs33941127"))

   mtx.test <- mtx.geno[1000:1020, 1:15]
   checkEquals(dim(mtx.test), c(21, 15))
   mtx <- skatifyGenotypeMatrix(mtx.test)
   checkEquals(dim(mtx), dim(mtx.test))
   checkEquals(as.numeric(table(mtx)), c(180, 100, 35))

      # now the entire mtx.geno
   mtx <- skatifyGenotypeMatrix(mtx.geno)
   checkEquals(dim(mtx), dim(mtx.geno))
                                      #     0    1   2
   checkEquals(as.numeric(table(mtx)), c(15251, 5864, 754))

   checkEquals(head(colnames(mtx), n=3), head(colnames(mtx.geno), n=3))
   checkEquals(head(rownames(mtx), n=3), head(rownames(mtx.geno), n=3))

} # test_skatifyGenotypeMatrix
#------------------------------------------------------------------------------------------------------------------------
ndufs2 <- function()
{
   tbl.cov <- createCovariatesTable()
   mtx.geno <- createGeno(pval.cutoff=0.05)
   save(mtx.geno, file=sprintf("mtx.geno.%dx%d.%s.RData", nrow(mtx.geno), ncol(mtx.geno),
                               gsub(" ", ".", Sys.time(), fixed=TRUE)))
   mtx.geno <- get(load("mtx.geno.1151x336.2021-10-04.08:38:50.RData"))
   dim(mtx.geno)

   shared.ids <- intersect(rownames(mtx.geno), tbl.cov$individualID)
   length(shared.ids) # 462
   mtx.geno <- mtx.geno[shared.ids,]
   tbl.cov <- subset(tbl.cov, individualID %in% shared.ids)
   tbl.cov$age_death[tbl.cov$age_death=="90+"] <- "90"
   tbl.cov$age_death <- as.numeric(tbl.cov$age_death)
      # 1143 without dx filtering, just 462 if dx 0 or 1 is used
   checkEquals(nrow(tbl.cov), 462)
   checkEquals(nrow(mtx.geno), 462)

     # gin up a strong signal, all patients with the first snp
   #has.rs429358 <- rep(0, nrow(tbl.cov))
   #has.rs429358[which(mtx.geno[, "19:45411941_T/C"] != "0/0")] <- 1
   #checkEquals(as.numeric(table(has.rs429358)), c(847, 296))
   #tbl.cov$has.rs429358 <- has.rs429358

   null.model <- SKAT_Null_Model(cogdx ~ v34 + v44 + apoe_genotype + msex  + age_death, tbl.cov)
   null.model <- SKAT_Null_Model(ceradsc ~ apoe_genotype, data=tbl.cov)
   null.model <- SKAT_Null_Model(cogdx ~ apoe_genotype + msex  + age_death, tbl.cov)
   null.model <- SKAT_Null_Model(ad ~ apoe_genotype + age_death + msex, tbl.cov)
   null.model <- SKAT_Null_Model(dx ~ apoe_genotype + age_death + msex, tbl.cov)
   #null.model <- SKAT_Null_Model(dx ~ apoe_genotype, tbl.cov)

   summary(lm(cogdx ~ apoe_genotype + msex , tbl.cov))
   summary(lm(cogdx ~ apoe_genotype, tbl.cov))
   summary(lm(ad ~ apoe_genotype + msex + age_death, tbl.cov))
   summary(lm(dx ~ apoe_genotype + msex + age_death, tbl.cov))$r.squared   # 0.03
   summary(lm(dx ~ apoe_genotype, tbl.cov))$r.squared                      # 0.017
   #summary(lm(dx ~ apoe_genotype + , tbl.cov))$r.squared                      # 0.017

   mtx.skat <- skatifyGenotypeMatrix(mtx.geno)
   x <- SKAT(mtx.skat, null.model)
   x$p.value
   x.0 <- SKAT(mtx.skat, null.model, method="SKATO")
   x.0$p.value

   rsids <- colnames(mtx.geno)
   tbl.rsid.locs <- subset(tbl.eqtl, rsid %in% rsids & study=="ampad-rosmap" &
                                     genesymbol=="NDUFS2")[, c("chrom", "hg38", "rsid")]
   tbl.rsid.locs <- tbl.rsid.locs[order(tbl.rsid.locs$hg38, decreasing=FALSE),]
   deleters <- which(tbl.rsid.locs$hg38 < 1)
   if(length(deleters) > 0)
       tbl.rsid.locs <- tbl.rsid.locs[-deleters,]

   with(tbl.rsid.locs, max(hg38) - min(hg38))     # 162186500
     # choose 4 snps at a time, overlapping by 3
   tbls.skat <- list()
   last.start <- ncol(mtx.skat) - 3
   sets <- seq_len(last.start)
   for(i in sets){
       range <- i:(i+3)
       printf("%d) range: %s", i, paste(range, collapse=" "))
       printf("   subsetting tbl.rsid.locs by range")
       rsids <- tbl.rsid.locs$rsid[range]
       rsids.valid <- intersect(rsids, colnames(mtx.skat))
       #printf("   subsetting mtx.skat by %d rsids.valid", length(rsids.valid))
       if(length(rsids.valid) == 0) next;
       mtx.skat.i <- mtx.skat[,rsids.valid, drop=FALSE]
       #printf("   mtx.skat.i:   %d x %d", nrow(mtx.skat.i), ncol(mtx.skat.i))
       x <- SKAT(mtx.skat.i, null.model, method="SKATO")
       if(x$p.value < 0.1)
           printf("%d: %f", i, x$p.value)
       score <- -log10(x$p.value)
       tbls.skat[[i]] <- data.frame(rsid=paste(rsids, collapse=","),
                                   score = score,
                                   pval = x$p.value,
                                   stringsAsFactors=FALSE)
       } # for i
   tbl.skat <- do.call(rbind, tbls.skat)
   tbl.skat <- tbl.skat[order(tbl.skat$score, decreasing=TRUE),]
   tbl.freq <- as.data.frame(sort(table(unlist(strsplit(head(tbl.skat, n=30)$rsid, split=","))),
                                  decreasing=TRUE))
   subset(tbl.freq, Freq >= 3)

      #-------------------------------------------
      # now with random groups of 4
      #-------------------------------------------

   tbl.skat <- list()
   for(i in 1:10000){
       range <- sort(sample(seq_len(ncol(mtx.skat)), size=4))
       printf("%d) range: %s", i, paste(range, collapse=" "))
       printf("   subsetting tbl.rsid.locs by range")
       #rsids <- tbl.rsid.locs$rsid[range]
       rsids <- colnames(mtx.skat)[range]
       if(any(is.na(rsids))) browser()
       rsids.valid <- intersect(rsids, colnames(mtx.skat))
       #printf("   subsetting mtx.skat by %d rsids.valid", length(rsids.valid))
       if(length(rsids.valid) == 0) next;
       mtx.skat.i <- mtx.skat[,rsids.valid, drop=FALSE]
       #printf("   mtx.skat.i:   %d x %d", nrow(mtx.skat.i), ncol(mtx.skat.i))
       x <- SKAT(mtx.skat.i, null.model, method="SKATO")
       if(x$p.value < 0.1)
           printf("%d: %s, %f", i, paste(rsids.valid, collapse=","), x$p.value)
       score <- -log10(x$p.value)
       tbls.skat[[i]] <- data.frame(rsid=paste(rsids, collapse=","),
                                   score = score,
                                   pval = x$p.value,
                                   stringsAsFactors=FALSE)
       } # for i
   tbl.skat <- do.call(rbind, tbls.skat)
   tbl.skat <- tbl.skat[order(tbl.skat$score, decreasing=TRUE),]
   tbl.freq <- as.data.frame(sort(table(unlist(strsplit(head(tbl.skat, n=30)$rsid, split=","))),
                                  decreasing=TRUE))
   subset(tbl.freq, Freq >= 3)
   tbl.rich <- subset(tbl.skat, pval <= 0.08)
   dim(tbl.rich)
   tbl.freq <- as.data.frame(sort(table(unlist(strsplit(tbl.rich$rsid, split=","))),
                                  decreasing=TRUE))
  head(tbl.freq, n=30)
   rsids.oi <- c("rs4575098", as.character(head(tbl.freq$Var1, n=3)))
   tbl.rich.eqtls <- subset(tbl.eqtl, rsid %in% rsids.oi)

   tbls.skat <- list()
   for(i in 1:ncol(mtx.skat)){
       mtx.skat.i <- mtx.skat[,i,drop=FALSE]
       x <- SKAT(mtx.skat.i, null.model, method="SKATO")
       #if(x$p.value < 0.1)
       #    printf("%d: %f", i, x$p.value)
       score <- -log10(x$p.value)
       tbls.skat[[i]] <- data.frame(rsid = colnames(mtx.skat)[i],
                                   score = score,
                                   pval = x$p.value,
                                   stringsAsFactors=FALSE)
       } # for i

   tbl.skat <- do.call(rbind, tbls.skat)
   tbl.skat <- tbl.skat[order(tbl.skat$score, decreasing=TRUE),]
   rownames(tbl.skat) <- NULL
   head(tbl.skat, n=10)
   tbl.logicFS <- run.logicFS(mtx.geno, tbl.cov, tbl.skat$rsid[1:10])

   etx <- EndophenotypeExplorer$new("NDUFS2", "hg38")
   mtx.rna <- etx$get.rna.matrix("old-rosmap")
   mtx.rna <- etx$get.rna.matrix("sage-eqtl-rosmap")
   mtx.rna[1:10, 1:10]
   x1 <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.rna, "rs34731114", study.name="rosmap")
   names(x1)
   x1$genotypes.rna
   x2$genotypes.vcf
   t.test(x1$wt["NDUFS2",], x1$mut["NDUFS2",])
   boxplot(x1$wt["NDUFS2",], x1$mut["NDUFS2",])

   x2 <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.rna, "rs3813623", study.name="rosmap")
   names(x2)
   x2$genotypes.rna
   x2$genotypes.vcf
   t.test(x2$wt["NDUFS2",], x2$mut["NDUFS2",])
   box2plot(x2$wt["NDUFS2",], x2$mut["NDUFS2",])

   subset(tbl.eqtl, rsid %in% c("rs4575098", "rs34731114", "rs3813623") & study=="ampad-rosmap" & pvalue < 1e-5)
    #   chrom      hg19      hg38       rsid      pvalue            ensg genesymbol        study tissue   assay
    #    chr1 161210576 161240786 rs34731114 7.91659e-06 ENSG00000158864     NDUFS2 ampad-rosmap  dlpfc unknown
    #    chr1 161169165 161199375  rs3813623 1.00585e-07 ENSG00000158864     NDUFS2 ampad-rosmap  dlpfc unknown
   subset(tbl.skat, rsid %in% c("rs34731114", "rs3813623"))
    # rank       rsid     score      pval
    #    3  rs3813623 1.3423514 0.0454620
    #    9 rs34731114 0.6970703 0.2008768
   etx$getAggregatedAlleleFrequencies("rs3813623")
    #      rsid ref       population A    T     G total A.freq    T.freq   G.freq  min.freq
    # rs3813623   G         European 0 4909 36315 41224      0 11.908112 88.09189 11.908112
    # rs3813623   G   African Others 0    5   183   188      0  2.659574 97.34043  2.659574
    # rs3813623   G       East Asian 0   93   251   344      0 27.034884 72.96512 27.034884

   tbls.freq <- lapply(tbl.skat$rsid[1:10], function(rsid) etx$getAggregatedAlleleFrequencies(rsid, quiet=FALSE)[1,])
   tbl.freq <- do.call(rbind.fill, tbls.freq)
   tbl.freq

   tbl.summary <- merge(tbl.skat, tbl.freq, by="rsid")

      # with    null.model <- SKAT_Null_Model(dx ~ apoe_genotype, tbl.cov)
   top.hits <- head(tbl.skat[order(tbl.skat$score, decreasing=TRUE),], n=10)$rsid
   top.hits
   eqtl.threshold <- 1e-4
   table(subset(tbl.eqtl, rsid %in% top.hits & study=="ampad-rosmap" & pvalue <  eqtl.threshold)$genesymbol)

   print(load("~/github/TrenaProjectAD/explore/rs4575098/ndufs2/motif.breaks.2119.variantsno.pvals.Sun-Sep-26-10:10:16-2021.RData"))
      # tbl.breaks, results
   tbl.breaks$pctDelta <- tbl.breaks$pctAlt - tbl.breaks$pctRef
   tbl.1127 <- subset(tbl.breaks, SNP_id == "rs33941127")
   dim(tbl.1127)  # 271 25
   #tbl.1127$pctDelta <- tbl.1127$pctAlt - tbl.1127$pctRef

   subset(tbl.1127, abs(pctDelta) > 0.15 & (pctAlt > 0.8 | pctRef > 0.8))[, c("geneSymbol", "pctAlt", "pctRef", "pctDelta")]
     # geneSymbol    pctAlt    pctRef   pctDelta
     #       EGR1 0.6242863 0.8059675 -0.1816812
     #        FOS 0.8310875 0.6713430  0.1597445
     #   FOS::JUN 0.8394816 0.6467574  0.1927243
     #       FOSB 0.8060674 0.6282966  0.1777708
     #      FOSL1 0.8360297 0.6853538  0.1506759
     #      FOSL1 0.8135891 0.6601670  0.1534220
     #      FOSL2 0.8197272 0.6693133  0.1504139
     #       JUNB 0.8047512 0.6464127  0.1583385
     #       JUNB 0.8554310 0.6945577  0.1608733
     #       JUND 0.8445769 0.6829076  0.1616693
     #     NKX2-3 0.8140993 0.6391464  0.1749528
     #     NKX3-2 0.8055512 0.6024556  0.2030956

} # dufs2
#----------------------------------------------------------------------------------------------------
run.logicFS <- function(mtx.geno, tbl.cov, rsids)
{
   mtx.char <- mtx.geno[, rsids]
   mtx <- matrix(0, nrow=nrow(mtx.char), ncol=ncol(mtx.char),
                 dimnames=list(rownames(mtx.char), colnames(mtx.char)))
   mtx[mtx.char == "0/1"] <- 1
   mtx[mtx.char == "1/1"] <- 1
   indices <- match(rownames(mtx), rownames(tbl.cov))
   labels <- tbl.cov$dx[indices]
   names(labels) <- rownames(tbl.cov)[indices]
   logicFS(mtx, labels,  B = 20, nleaves = 10, rand = 1234, anneal.control = my.anneal)

} # logicFS
#----------------------------------------------------------------------------------------------------
skat.with.all.combinations.of.topRanked.eqtls <- function()
{
   rsids <- subset(tbl.eqtl, pvalue <= 0.00001)$rsid
   length(rsids)
   rsids <- intersect(rsids, colnames(mtx.skat))
   length(rsids)
   rsid.sets <- combn(rsids, 4, simplify=FALSE)
   length(rsid.sets)

   tbls.skat <- list()
   i <- 0
   for(set in rsid.sets){
      i <- i + 1
      mtx.skat.i <- mtx.skat[,set]
      x <- SKAT(mtx.skat.i, null.model, method="SKATO")
      if(x$p.value < 0.1)
          printf("%d: %s, %f", i, paste(set, collapse=","), x$p.value)
      score <- -log10(x$p.value)
      tbls.skat[[i]] <- data.frame(rsid=paste(rsids, collapse=","),
                                   score = score,
                                   pval = x$p.value,
                                   stringsAsFactors=FALSE)
       } # for i
   tbl.skat <- do.call(rbind, tbls.skat)
   tbl.skat <- tbl.skat[order(tbl.skat$score, decreasing=TRUE),]
   tbl.freq <- as.data.frame(sort(table(unlist(strsplit(head(tbl.skat, n=30)$rsid, split=","))),
                                  decreasing=TRUE))


} # skat.with.all.combinations.of.topRanked.eqtls
#----------------------------------------------------------------------------------------------------

