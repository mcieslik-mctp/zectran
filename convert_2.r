hgnc2ucsc = function(hgnc) {
    out = suppressWarnings(
        AnnotationDbi::select(org.Hs.eg.db, keys=c(hgnc), columns=c("UCSCKG"), keytype="SYMBOL")
        )
    colnames(out) = c("hgnc", "ucsc")
    return(out)
}

ucsc2grl = function(ucsc) {
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
        protein_id=NA
        )
    tx_gfull = suppressWarnings(
        keepSeqlevels(tx_gfull, CCHR)
        )
    if (length(tx_gfull) > 0) {
        seqlengths(tx_gfull) = CHR_LEN[names(seqlengths(tx_gfull))]
        tx_grs = split(tx_gfull, tx_gfull$transcript_id)
    } else {
        tx_grs = tx_gfull # empty range
    }
    return(tx_grs)
}

grl2probes = function {
    findSpanOverlaps(grl, GR450K)
}
