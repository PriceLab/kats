library(VariantAnnotation)
# apoe:  chr19:45,408,053-45,413,650
big.region <- GRanges(seqnames="19", IRanges(start=45360000, end=45460000))
range(big.region)
width(big.region)
vcf.big <- readVcf(vcfFile, "hg19", big.region)
length(vcf.big)  # 2754
file.name <- "chr19-apoe-100k-hg19.vcf"
gz.filename <- sprintf("%s.gz", file.name)
writeVcf(vcf.big, filename=file.name)

bgzip(file.name, dest=gz.filename, overwrite=TRUE)
indexTabix(gz.filename, format="vcf")
grep(file.name, list.files(), value=TRUE)
  # "chr19-apoe-100k-hg19.vcf"        "chr19-apoe-100k-hg19.vcf.gz"     "chr19-apoe-100k-hg19.vcf.gz.tbi"
