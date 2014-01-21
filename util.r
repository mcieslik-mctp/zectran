NULL_GR = GRanges()
NULL_GRL = GRangesList()

stringToRange = function(string) {
    tryCatch({
        tmp = str_replace_all(string, ",", "")
        fs = unlist(str_split(str_split(tmp, pattern=":")[[1]], "-"))
        GRanges(seqnames=fs[1], ranges=IRanges(start=as.integer(fs[2]), end=as.integer(fs[3])), strand="*", type="range")
    }, error=function(e) {NULL_GR})
}

findSpanOverlaps = function(query, subject) {
    if (length(query) > 0) {
        gr = unlist(query)
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

getEmptyRange = function(transcript) {
    NULL_GR
}

getEmptyTranscript = function(transcripts) {
    NULL_GRL
}

getFirstTranscript = function(transcripts) {
    if (length(transcripts) == 0) {
        return(getEmptyTranscript())
    }
    transcripts[[1]]
}

getUnionTranscript = function(transcripts) {
    if (length(transcripts) == 0) {
        return(getEmptyTranscript())
    }
    all = unlist(transcripts)
    uni = reduce(all)
    mcols(uni)$type = "exon"
    mcols(uni)$exon_number = ifelse(strand(uni) == "+", 1:length(uni), length(uni):1)
    unis = GRangesList()
    gene_id = unique(mcols(all)$gene_id)
    mcols(uni)$gene_id = gene_id
    unis[[gene_id]] = uni
    return(unis)
}

getCanonicalTranscript = function(transcripts) {
    ## transcript with the longest exons
    if (length(transcripts) == 0) {
        return(getEmptyTranscript())
    }
    gene_id = unique(mcols(transcripts)$gene_id)
    idx = which.max(sapply(width(transcripts), sum))
    sel = transcripts[idx:idx]
    return(sel)
}

flankTranscript = function(transcript, upstream, downstream) {
    if (length(transcript) == 0) {
        return(getEmptyRange())
    }
    print(transcript)
    exon_number = mcols(transcript)$exon_number
    print(exon_number)
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

silentExec = function(...) {
    capture.output(
        {res = suppressMessages(suppressWarnings(...))})
    return(res)
}

safeEval = function(f, ...) {
    tryCatch(do.call(f, list(...)), error=function(e) NULL)
}

safeSwitch = function(...) {
    tryCatch(switch(...), error=function(e) NULL)
}

