library(podkat)
library(RUnit)
#----------------------------------------------------------------------------------------------------
file.dir <- "~/github/TrenaProjectAD/explore/WGS-Harmonization/metadata-all"

specimen.file <- file.path(file.dir, "ROSMAP_biospecimen_metadata.csv")
checkTrue(file.exists(specimen.file))
tbl.spec <- read.table(specimen.file, sep=",", header=TRUE, nrow=-1)

clinical.file <- file.path(file.dir, "ROSMAP_clinical.csv")
checkTrue(file.exists(clinical.file))
tbl.clin <- read.table(clinical.file, sep=",", header=TRUE, nrow=-1)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{

} # runTests
#----------------------------------------------------------------------------------------------------
test_createCovariatesTable <- function()
{

   tbl.cov <- subset(tbl.clin, Study=="ROS")[, c("individualID", "apoe_genotype", "cogdx", "age_death", "msex")]
   dim(tbl.cov)  # 1451 6
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
    # reduce these to 0 1 1 5 4 2

   tbl.cov$age_death[tbl.cov$age_death=="90+"] <- "101"

   tbl.cov$age_death <- round(as.numeric(tbl.cov$age_death), digits=0)
   outcomes <- tbl.cov$cogdx[keepers]
   apoe <- tbl.cov$apoe_genotype[keepers]
   age <- tbl.cov$age_death[keepers]

   null.model <- SKAT_Null_Model(outcomes ~ apoe + age, out_type="C")



} # test_createCovariatesTable
#----------------------------------------------------------------------------------------------------
