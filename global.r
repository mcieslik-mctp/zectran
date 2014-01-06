options(stringsAsFactors=FALSE)
SYNC = Sys.getenv("SYNC")

## imports
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg19")) # Hsapiens
suppressWarnings(library(org.Hs.eg.db)) # conversions
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressMessages(library("GenomicFeatures")) # GenomicRanges

## CHROMOSOMES
ICHR = c(1:22, "X", "Y")
CCHR = paste("chr", ICHR, sep="")
CHR_LEN = seqlengths(Hsapiens)[CCHR]

## TXDB
get_TX = function(model="ensembl") {
    TXDB = switch(model,
        "ensembl"={
            TXS_FN = file.path(Sys.getenv("SYNC"), "data",
                "txdb-ensembl_73_GRCh37.p12.sqlite")
            loadDb(TXS_FN)
        },
        "ucsc"={
            TXS_FN = file.path(Sys.getenv("SYNC"), "data",
                "txdb-ucsc_19-09-2013_hg19.sqlite")
            loadDb(TXS_FN)
        })
}

TX_ENSEMBL_73 = get_TX("ensembl")
TX_UCSC_HG19 = get_TX("ucsc")

## 450K
GR_450K = sort(readRDS("tables/FDb.InfiniumMethylation.hg19_09-10-2013.rds"))
GR_450K_IDS = names(GR_450K)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ANN_450K = IlluminaHumanMethylation450kanno.ilmn12.hg19@data
