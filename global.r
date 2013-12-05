##
SYNC = Sys.getenv("SYNC")

source(file.path(SYNC, "libs", "lib.r"))
source("hdf5.r")
source("convert.r")

library(h5r)
library(stringr)
library(shiny)
library(ggbio) ## ggplot
library(biomaRt)
suppressMessages(library(data.table))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg19")) # Hsapiens
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressMessages(library("GenomicFeatures")) # GenomicRanges 

## CHROMOSOMES
ICHR = c(1:22, "X", "Y")
CCHR = paste("chr", ICHR, sep="")
CHR_LEN = seqlengths(Hsapiens)[CCHR] # hg19 lengths

##
HGNC2ENSG = data.table(fread("tables/ensg2hgnc.tbl"))
setnames(HGNC2ENSG,colnames(HGNC2ENSG),c("ensg", "hgnc"))
setkey(HGNC2ENSG, "hgnc")

## TXDB
get_TXDB = function() {
    TXS_FN = file.path(Sys.getenv("SYNC"), "data", "txdb-ensembl_73_GRCh37.p12.sqlite")
    TXDB = loadDb(TXS_FN)
    return(TXDB)
}
TXDB = get_TXDB()

## H5FILE
H5DATA = "tables/450k_tbl.h5"
H5FILE = H5File(H5DATA, "r")
H5DATA_ROW = "tables/450k_tbl_row.h5"
H5FILE_ROW = H5File(H5DATA_ROW, "r")

## TCGA META
ALIQUOT_FN = "tables/bio:bcr_aliquot_uuid.tsv"
SDRF_FN = "tables/mage-450k_merged_sdrf.txt"

make_aliquot_meta = function(aliquot_fn=ALIQUOT_FN, sdrf_fn=SDRF_FN) {
    aliquot = tread(aliquot_fn)
    aliquot = aliquot[!is.na(aliquot$bio.bcr_aliquot_uuid),]
    ALIQUOT = aliquot[,c("bio.bcr_aliquot_uuid", "bio.bcr_sample_uuid", "admin.disease_code", "bio.sample_type_id")]
    ALIQUOT$bio.bcr_aliquot_uuid = toupper(ALIQUOT$bio.bcr_aliquot_uuid)
    ALIQUOT$bio.bcr_sample_uuid = toupper(ALIQUOT$bio.bcr_sample_uuid)
    
    sdrf = tread(sdrf_fn)
    ARRAY = sdrf[c("Extract.Name", "Array.Data.File")]
    names(ARRAY) = c("bio.bcr_aliquot_uuid", "file")
    ARRAY$bio.bcr_aliquot_uuid = toupper(ARRAY$bio.bcr_aliquot_uuid)

    ARRAY$array_450k_id = str_sub(ARRAY$file, end=-10)
    ARRAY$file = NULL
    ARRAY = unique(ARRAY)
    ALIMETA = merge(ALIQUOT, ARRAY, by="bio.bcr_aliquot_uuid")
}
ALIMETA = make_aliquot_meta()
TSAMPLES = ALIMETA[ALIMETA$bio.sample_type_id <  10, "array_450k_id"]
NSAMPLES = ALIMETA[ALIMETA$bio.sample_type_id >= 10, "array_450k_id"]

## 
GR450K = sort(readRDS("tables/FDb.InfiniumMethylation.hg19_09-10-2013.rds"))
GR450K_IDS = names(GR450K)

## minfi
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ANN450K = IlluminaHumanMethylation450kanno.ilmn12.hg19@data

COHORT_DESC = tread("tables/diseaseStudy.txt")
rownames(COHORT_DESC) = COHORT_DESC$Study.Abbreviation
colnames(COHORT_DESC) = c("short", "full")

print("global")
## Genome Annotation
## ENSEMBL_GTF = readRDS("tables/Homo_sapiens.GRCh37.73.rds")
## ENSEMBL_GTF = keepSeqlevels(ENSEMBL_GTF, c(1:22, "X", "Y"))
## ENSEMBL_GTF = renameSeqlevels(ENSEMBL_GTF, paste("chr", seqlevels(ENSEMBL_GTF), sep=""))

## HUGO2ENSG = data.table(as.data.frame(unique(mcols(ENSEMBL_GTF)[,c("gene_name", "gene_id")])))
## setkey(HUGO2ENSG, gene_name)
## ENSG2ENST = data.table(as.data.frame(unique(mcols(ENSEMBL_GTF)[,c("gene_id", "transcript_id")])))
## setkey(ENSG2ENST, gene_id)
