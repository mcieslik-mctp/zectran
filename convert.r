##
hgnc2ensg = function(hgnc) {
    HGNC2ENSG[hgnc]
}

## print("2ensg-1")
## out = getBM(MART, attributes=c("hgnc_symbol", "ensembl_gene_id"), filters=c("hgnc_symbol"), values=list(hgnc))
## print("2ensg-2")
## colnames(out) = c("hgnc", "ensg")
## return(out)
## ## MART
## get_MART = function() {
##     marts = listMarts()
##     ensembl_version = str_extract(as.character(marts[marts$biomart == "ensembl","version"]), "([0-9]+)")
##     stopifnot(ensembl_version == "73")
##     MART = useMart("ensembl", "hsapiens_gene_ensembl")
##     return(MART)
## }
## MART = get_MART()
## out = getBM(MART, attributes=c("hgnc_symbol", "ensembl_gene_id"))


##
ensg2enst = function(ensg) {
    out = suppressWarnings(
        AnnotationDbi::select(TXDB, columns=c("GENEID", "TXNAME"), keytype="GENEID", keys=ensg))
    colnames(out) = c("ensg", "enst")
    return(out)
}

##
enst2grl = function(enst) {
    tx_query = suppressWarnings(
        AnnotationDbi::select(TXDB, columns=columns(TXDB), keytype="TXNAME", keys=enst))
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
        keepSeqlevels(tx_gfull, c(1:22, "X", "Y"))
        )
    if (length(tx_gfull) > 0) {
        tx_gfull = renameSeqlevels(tx_gfull, paste("chr", seqlevels(tx_gfull), sep=""))
        seqlengths(tx_gfull) = CHR_LEN[names(seqlengths(tx_gfull))]
        tx_grs = split(tx_gfull, tx_gfull$transcript_id)
    } else {
        tx_grs = tx_gfull # empty range
    }
    return(tx_grs)
}

grl2probes = function(grl) {
    gr = unlist(grl)
    if (length(gr) > 0) {
        query = GRanges(
        seqnames(gr)[1],
        IRanges(min(start(gr)), max(end(gr))),
        strand="*"
        )
        probes = GR450K[findOverlaps(query, GR450K)@subjectHits]
    } else {
        probes = GR450K[NULL]
    }
    return(probes)
}

gene_flank = function(transcript, upstream, downstream) {
    exon_number = mcols(transcript)$exon_number
    first_exon = transcript[exon_number == min(exon_number)]
    last_exon  = transcript[exon_number == max(exon_number)]
    if (as.vector(strand(first_exon)) == "+") {
        tss = start(first_exon) - upstream
        tts = end(last_exon) + downstream
        tss_tts = IRanges(start=tss, end=tts)
    } else {
        tss = end(first_exon) + upstream
        tts = start(last_exon) - downstream
        tss_tts = IRanges(start=tts, end=tss)
    }
    ##
    gr = GRanges(
            seqnames=seqnames(first_exon),
            ranges=tss_tts,
            strand=strand(first_exon),
            type="range",
            gene_id=mcols(first_exon)$gene_id
            #transcript_id=mcols(first_exon)$transcript_id[1]
            )
    return(gr)
}

print("convert")
