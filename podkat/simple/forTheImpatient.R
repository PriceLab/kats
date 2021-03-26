library(podkat)
library(VariantAnnotation)
phenoFileLin <- system.file("examples/example1lin.csv", package="podkat")
phenoFileLog <- system.file("examples/example1log.csv", package="podkat")
vcfFile <- system.file("examples/example1.vcf.gz", package="podkat")

pheno.c <- read.table(phenoFileLin, header=TRUE, sep=",")
dim(pheno.c)  # 200 3
head(pheno.c)
model.c <- nullModel(y ~ ., pheno.c)
class(model.c)
model.c

data(hgA)
class(hgA)  # GRanges
length(hgA) # 1
hgA         # invented genome, chr1 1-200000

windows <- partitionRegions(hgA, width=5000, overlap=0.5)

vcf <- readVcf(vcfFile)
mtx.geno <- geno(vcf)$GT  # [1] 3845  200

  # class: CollapsedVCF 
  # dim: 3845 200 
  # rowRanges(vcf):
  #   GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
  # info(vcf):
  #   DataFrame with 0 columns: 
  # geno(vcf):
  #   List of length 1: GT
  # geno(header(vcf)):
  #       Number Type   Description
  #    GT 1      String Genotype   

  # mtx.geno[1:10, 1:10]
  #        S1    S2    S3    S4    S5    S6    S7    S8    S9    S10  
  # snv:1  "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0"
  # snv:2  "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0"
  # snv:3  "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0"
  # snv:4  "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0"
  # snv:5  "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0"
  # snv:6  "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "1/0" "0/0" "0/0"
  # snv:7  "0/0" "0/0" "0/0" "0/0" "0/0" "1/0" "0/0" "0/0" "0/0" "0/0"
  # snv:8  "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0"
  # snv:9  "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0"
  # snv:10 "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0"



  #      readGenotypeMatrix(file, regions, subset,
  #                         noIndels=TRUE, onlyPass=TRUE,
  #                         na.limit=1, MAF.limit=1,
  #                         na.action=c("impute.major", "omit", "fail"),
  #                         MAF.action=c("invert", "omit", "ignore", "fail"),
  #                         sex=NULL)

geno <- readGenotypeMatrix(vcfFile)
dim(geno)  # 200 962

  # as.matrix(geno)[1:10, 1:10]
  #     snv:6 snv:7 snv:9 snv:11 snv:12 snv:18 snv:21 snv:23 snv:29 snv:30
  # S1      0     0     0      0      0      0      0      0      0      0
  # S2      0     0     0      0      0      0      1      0      0      0
  # S3      0     0     0      0      0      0      0      0      0      0
  # S4      0     0     0      0      0      1      0      0      0      0
  # S5      0     0     0      0      0      1      0      0      0      0
  # S6      0     1     0      0      0      0      0      0      0      1
  # S7      0     0     0      0      0      0      0      0      0      0
  # S8      1     0     0      0      0      0      0      0      0      0
  # S9      0     0     0      0      0      2      0      0      0      0
  # S10     0     0     0      0      0      1      0      0      0      0

res.c <- assocTest(geno, model.c, windows)
print(res.c)
plot(res.c, which="p.value")
  # Overview of association test:
  # 	Null model: linear 
  # 	Number of samples: 200 
  # 	Number of regions: 79 
  # 	Number of regions without variants: 0 
  # 	Average number of variants in regions: 24.1 
  # 	Genome: hgA 
  # 	Kernel: linear.podkat 
  # 	p-value adjustment: none 
  # 
  # Overview of significance of results:
  # 	Number of tests with p < 0.05: 8
  # 
  # Results for the 8 most significant regions:
  #   seqnames  start    end width  n         Q      p.value
  # 1     chr1   7501  12500  5000 31 769748.34 1.294084e-07
  # 2     chr1  10001  15000  5000 33 764828.81 4.874460e-06
  # 3     chr1 140001 145000  5000 15  79937.68 3.599077e-03
  # 4     chr1   5001  10000  5000 34 152555.30 9.785569e-03
  # 5     chr1 132501 137500  5000 21  89287.55 1.349559e-02
  # 6     chr1 142501 147500  5000 23  94629.68 3.338620e-02
  # 7     chr1  42501  47500  5000 19  58191.23 3.341032e-02
  # 8     chr1  25001  30000  5000 23 103713.12 3.754557e-02

windows <- partitionRegions(hgA, width=250, overlap=0.5)
res.c <- assocTest(geno, model.c, windows)
print(res.c)

res.c <- p.adjust(res.c)
print(res.c)
plot(res.c, which="p.value.adj")


   #------------------------------------------------------------
   # now the binary model
   #------------------------------------------------------------

pheno.b <- read.table(phenoFileLog, header=TRUE, sep=",")
model.b <- nullModel(y ~ ., pheno.b)
model.b

res.b <- assocTest(vcfFile, model.b, windows)
print(res.b)
res.b <- p.adjust(res.b)
print(res.b)
plot(res.b, which="p.value.adj")
