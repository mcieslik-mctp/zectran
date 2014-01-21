safeSelect = function(db, keys, columns, keytype, names) {
    out = tryCatch(
        as.data.frame(
            suppressWarnings(
                AnnotationDbi::select(db, keys=keys, columns=columns, keytype=keytype)
                )),
        error=function(e) data.frame(replicate(length(columns) + 1, character(0))))
    colnames(out) = names
    return(out)
}


hgnc2ucsc = function(hgnc) {
    safeSelect(org.Hs.eg.db, hgnc, columns=c("UCSCKG"), keytype="SYMBOL", c("hgnc", "ucsc"))
}

hgnc2name = function(hgnc) {
    safeSelect(org.Hs.eg.db, hgnc, columns=c("GENENAME"), keytype="SYMBOL", c("hgnc", "name"))
}

hgnc2eg = function(hgnc) {
    safeSelect(org.Hs.eg.db, hgnc, columns=c("ENTREZID"), keytype="SYMBOL", c("hgnc", "eg"))
}

ucsc2ucsc = function(ucsc) {
    patterns = sapply(str_split(ucsc, "\\."), "[", 1)
    unlist(lapply(patterns, function(pattern) {
        keys(TX_UCSC_HG19, keytype="TXNAME", pattern=pattern)
    }))
}

ucsc2grl = function(ucsc) {
    tryCatch({
        tx_query = suppressWarnings(
            AnnotationDbi::select(TX_UCSC_HG19, columns=columns(TX_UCSC_HG19), keytype="TXNAME", keys=ucsc))
        tx_ifull = IRanges(start=tx_query$EXONSTART, end=tx_query$EXONEND)
        tx_gfull = GRanges(seqnames=tx_query$TXCHROM, ranges=tx_ifull,
            strand=tx_query$TXSTRAND,
            source=NA,
            type="exon",
            score=NA,
            phase=NA,
            gene_id=tx_query$GENEID,
            transcript_id=tx_query$TXNAME,
            exon_number=tx_query$EXONRANK,
            gene_name=NA,
            gene_biotype=NA,
            transcript_name=NA,
            exon_id=tx_query$EXONNAME,
            protein_id=NA)
        tx_gfull = suppressWarnings(
            keepSeqlevels(tx_gfull, CCHR))
        if (length(tx_gfull) > 0) {
            seqlengths(tx_gfull) = CHR_LEN[names(seqlengths(tx_gfull))]
            split(tx_gfull, tx_gfull$transcript_id)
        } else {
            NULL_GR
        }
    }, error=function(e) NULL_GR)
}

grl2probes = function(grl) {
    findSpanOverlaps(grl, GR_450K)
}
