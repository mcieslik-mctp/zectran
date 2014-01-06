findSpanOverlaps = function(query, subject) {
    gr = unlist(grl)
    if (length(gr) > 0) {
        q = GRanges(
        seqnames(gr)[1],
        IRanges(min(start(gr)), max(end(gr))),
        strand="*"
        )
        probes = subject[findOverlaps(q, subject)@subjectHits]
    } else {
        probes = subject[NULL]
    }
    return(probes)
}

flankTranscript = function(transcript, upstream, downstream) {
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
    gr = GRanges(
            seqnames=seqnames(first_exon),
            ranges=tss_tts,
            strand=strand(first_exon),
            type="range",
            gene_id=mcols(first_exon)$gene_id
            )
    seqlengths(gr) = CHR_LEN[names(seqlengths(gr))]
    return(gr)
}
