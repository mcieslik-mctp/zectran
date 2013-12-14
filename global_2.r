options(stringsAsFactors=FALSE)
SYNC = Sys.getenv("SYNC")

suppressMessages(library("GenomicFeatures"))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg19")) # Hsapiens

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
