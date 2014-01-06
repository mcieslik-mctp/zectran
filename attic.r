## Tables
.cohort_description = function(cohorts) {
    if (length(cohorts) == 0) {
        return(NULL)}
    
    paste(cohorts, COHORT_DESC[cohorts, "full"], sep=": ", collapse="\n")
}

.cohort_samples = function(cohorts) {
    if (length(cohorts) == 0) {
        return(NULL)}
    
    df = ldply(cohorts, function(cohort) {
        samples = get_samples(cohort)
        cbind(cohort=cohort, array_450k_id=samples)
    })
    df = merge(df, ALIMETA, by="array_450k_id")
    df$status = ifelse(df$bio.sample_type_id < 10, "Tumor", "Normal")
    df$status = paste(df$status, " (", df$bio.sample_type_id, ") ", sep="")
    df$admin.disease_code = NULL
    df$bio.sample_type_id = NULL
    return(df)
}

.cohort_summary = function(cohorts) {
    if (length(cohorts) == 0) {
        return(NULL)
    }
    
    tab = as.matrix(t(xtabs(data=ALIMETA, rep(1,nrow(ALIMETA)) ~ admin.disease_code + bio.sample_type_id)))
    cols = colnames(tab)
    rows = rownames(tab)
    storage.mode(tab) = "integer"
    rownames(tab) = paste(c("primary", "recurrent", "blood", "new_primary", "metastasis", "tissue_normal"),
                " (", rows, ")", sep="")
    tab = as.data.frame.matrix(tab)
    tab = tab[,colnames(tab) %in% cohorts, drop=FALSE]
    return(tab)
}

.ensts_grl = function(ensts) {
    if (length(ensts) == 0) {
        return(NULL)}
    
    ensts_grl = enst2grl(ensts)
    return(ensts_grl)
}

.ensts_grl_select = function(ensts_grl, select) {
    if (length(ensts_grl) == 0) {
        return(NULL)}
    
    switch(select,
           "union"={
               ensts_all = unlist(ensts_grl)
               uni = reduce(ensts_all)
               mcols(uni)$type = "exon"
               mcols(uni)$exon_number = ifelse(strand(uni) == "+", 1:length(uni), length(uni):1)
               unis = GRangesList()
               gene_id = unique(mcols(ensts_all)$gene_id)
               mcols(uni)$gene_id = gene_id
               unis[[gene_id]] = uni
               unis
           },
           "longest"={
               ensts_width = unlist(lapply(ensts_grl, function(transcript) {
                   max(end(transcript)) - min(start(transcript))
               }))
               select = names(ensts_width[order(ensts_width, decreasing=TRUE)[1]])
               ensts_grl[select]
           },
           {ensts_grl[select]
           })
}

.ensts_grl_select_region = function(ensts_grl_select, region) {
    if (length(ensts_grl_select) == 0) {
        return(NULL)
    }
    
    transcript = ensts_grl_select[[1]] # should be only one
    exon_number = mcols(transcript)$exon_number
    switch(region,
           "TSS slim (500, 200)"={
               first_exon = transcript[exon_number == 1]
               promoters(first_exon, 500, 200)
           },
           "TSS wide (2000, 500)"={
               first_exon = transcript[exon_number == 1]
               promoters(first_exon, 2000, 500)
           },
           "gene regulatory region"={
               gene_flank(transcript, 2000, 2000)
           }
           )
}

.ensts_grl_select_region_str = function(ensts_grl_select_region) {
    if (length(ensts_grl_select_region) == 0) {
        return(NULL)}
    
    sprintf("%s:%s-%s", seqnames(ensts_grl_select_region),
            start(ensts_grl_select_region),
            end(ensts_grl_select_region))
}

.ensts_grl_select_region_str_parse = function(ensts_grl_select_region_str) {
    if ((    length(ensts_grl_select_region_str) == 0 ) ||
        (str_length(ensts_grl_select_region_str) == 0 ) ) {
        return(NULL)}
    
    str_replace_all(ensts_grl_select_region_str, ",", "")
    fields = unlist(str_split(str_split(ensts_grl_select_region_str, pattern=":")[[1]], "-"))
    chr = fields[1]
    start = as.integer(fields[2])
    end = as.integer(fields[3])
    gr = GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), strand="*", type="range")
    return(gr)
}

.cpg_info = function(cpg) {
    if (length(cpg) == 0) {
        return(NULL)
    }
    
    rng = GR450K[cpg,]
    info = as(rng, "data.frame")
    info = info[,!(names(info) %in% "sourceSeq")]
    info = rename(info, c(seqnames="chr"))
    return(info)
}

.cpg_islands_info = function(cpg) {
    if (length(cpg) == 0) {
        return(NULL)}

    islands_info = as.data.frame(ANN450K$Islands.UCSC[cpg,])
    return(islands_info)
}

## Plots
.ensts_tracks_gq = function(ensts_grl) {
    if (length(ensts_grl) == 0) {
        return(NULL)}
    
    plt = tracks(
        autoplot(ensts_grl, group.selfish=TRUE, gap.geom="arrow") + theme_bw() + scale_x_sequnit("kb")
        )
    return(plt)
}

.ensts_select_tracks_gq = function(ensts_select_grl, ensts_select_grl_region, ensts_select_grl_region_probes) {
    if ((length(ensts_select_grl) == 0) ||
        (length(ensts_select_grl_region) == 0) ||
        (length(ensts_select_grl_region_probes) == 0)
        ) {return(NULL)}
    
    plt = tracks(
        autoplot(ensts_select_grl, group.selfish=TRUE, gap.geom="arrow") + theme_bw() + scale_x_sequnit("kb"),
        autoplot(ensts_select_grl_region, group.selfish=TRUE) + theme_bw() + scale_x_sequnit("kb"),
        autoplot(ensts_select_grl_region_probes, group.selfish=TRUE) + theme_bw() + scale_x_sequnit("kb"),
        heights=c(0.6,0.2,0.2)
        )
    return(plt)
}

.cg_boxplot_gq = function(probes_betas, comparison, title) {
    if (length(probes_betas) == 0) {
        return(NULL)}
    df = ldply(probes_betas, function(betas) {
        melt(betas)
    })
    names(df) = c("cohort", "cg", "sample", "beta")
    df$status = ifelse(df$sample %in% TSAMPLES, "tumor", "normal")
    df$cg = ordered(df$cg, rownames(probes_betas[[1]]))
    base_plt = ggplot() + aes(x=cg, y=beta) + theme_bw() +
            geom_boxplot() +
            scale_fill_discrete() +
            scale_color_discrete(guide=FALSE) +
            ggtitle(title) +
            theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
    if (comparison == "tumor") {
        plt = (base_plt + facet_grid(cohort ~ .)) %+% df[df$status == "tumor",]
    } else if (comparison == "normal") {
        plt = (base_plt + facet_grid(cohort ~ .)) %+% df[df$status == "normal",]
    } else if (comparison == "tumor-normal") {
        plt = (base_plt %+% df) + aes(fill=status, color=status) + facet_grid(cohort ~ .)  
    } else {
        plt = NULL
    }
    return(plt)
}

.cg_corrplot_gq = function(probes_betas, title) {
    if (length(probes_betas) == 0) {
        return(NULL)}
    
    beta_mx = data.frame(t(splat(cbind)(probes_betas)))
    beta_corr = melt(cor(beta_mx, use="complete"))
    beta_corr$X1 = ordered(beta_corr$X1, colnames(beta_mx))
    beta_corr$X2 = ordered(beta_corr$X2, colnames(beta_mx))
    plt = ggplot(beta_corr) + aes(x=X1, y=X2, fill=value) + theme_bw() +
        geom_tile() + 
        scale_fill_gradient2(high="red", mid="black", low="blue", guide=FALSE) + xlab(NULL) + ylab(NULL) +
        ggtitle(title) +
        theme(axis.ticks.x = element_blank(),
               axis.text.x = element_blank())
    return(plt)
}

.ensts_grl_select_region_limit_probes_names = function(probes) {
    if (length(probes) == 0) {
        return(list())} # ?
    
    names(probes)
}

.plot_pca_sa = function() {
    library(mvtnorm)

    x1 = runif(1)
    y1 = runif(1)
    xy1 <- rmvnorm(n=100, mean=c(rnorm(1), rnorm(1)), sigma=matrix(c(1, 0.2*y1, 0.2*y1, 1), nrow=2))
    x2 = runif(1)
    y2 = runif(1)
    xy2 <- rmvnorm(n=100, mean=c(rnorm(1), rnorm(1)), sigma=matrix(c(1, 0.2*y2, 0.2*y2, 1), nrow=2))
    df = data.frame(rbind(xy1, xy2))
    colnames(df) = c("pca1", "pca2")
    df$status = c(rep("T", 100), rep("N", 100))
    plt = ggplot(df) + geom_point() + theme_bw() + xlab("PCA 1") + ylab("PCA 2") + aes(x=pca1, y=pca2, color=status)
    return(plt)
    
    
}

ensts = reactive({
    if (is.null(input$select_ensg_gq) || input$select_ensg_gq == "invalid gene") {
        return(NULL)
    }
    ensg = input$select_ensg_gq
    ensts = ensg2enst(ensg)$enst
    return(ensts)
})

output$select_enst_gq = renderUI(
    selectInput("select_enst_gq", label="select transcript", choices=c("union", "longest", ensts()))
    )

ensts_grl = reactive({.ensts_grl(ensts())})

ensts_tracks_gq = reactive({
    input$plot_ensg_gq
    isolate({
        tks = .ensts_tracks_gq(ensts_grl())
        print(tks)
    })})
output$ensts_tracks_gq = renderPlot(ensts_tracks_gq())


ensts_grl_select = reactive({.ensts_grl_select(ensts_grl(), input$select_enst_gq)})

ensts_grl_select_region = reactive({.ensts_grl_select_region(ensts_grl_select(), input$region_select_gq)})

ensts_grl_select_region_str = reactive({.ensts_grl_select_region_str(ensts_grl_select_region())})

output$query_region_limit_gq = renderUI(
    textInput("query_region_limit_gq", label="limit region", ensts_grl_select_region_str()))

ensts_grl_select_region_limit = reactive({.ensts_grl_select_region_str_parse(input$query_region_limit_gq)})

ensts_grl_select_region_limit_probes = reactive({grl2probes(ensts_grl_select_region_limit())})

ensts_select_tracks_gq = reactive({
    input$plot_enst_gq
    isolate({
        tks = .ensts_select_tracks_gq(
            ensts_grl_select(),
            ensts_grl_select_region_limit(),
            ensts_grl_select_region_limit_probes())
        print(tks)
    })})
output$ensts_select_tracks_gq = renderPlot(ensts_select_tracks_gq())

ensts_grl_select_region_limit_probes_names = reactive({.ensts_grl_select_region_limit_probes_names(
    ensts_grl_select_region_limit_probes())})

output$select_cg_gq = renderUI(selectInput("select_cg_gq", label="select probes",
    choices=ensts_grl_select_region_limit_probes_names(), multiple=TRUE))

## Analyses
## cg info
cpg_info = reactive({
    .cpg_info(input$select_cg_gq)
})
output$cpg_info = renderTable(cpg_info())

cpg_islands_info = reactive({
    .cpg_islands_info(input$select_cg_gq)
})
output$cpg_islands_info = renderTable(cpg_islands_info())

## w/ beta values
probes_betas = reactive({
    get_beta_n(input$cohort_select_gq, input$select_cg_gq)
})
title = reactive({
    sprintf("%s (%s) [%s]", input$query_hgnc_gq, input$select_ensg_gq, ensts_grl_select_region_str())
})

cg_boxplot_gq = reactive({
    plt = .cg_boxplot_gq(probes_betas(), comparison=input$comparison_select_gq, title=title())
    print(plt)
})
output$cg_boxplot_gq = renderPlot(cg_boxplot_gq())

cg_corrplot_gq = reactive({
    plt = .cg_corrplot_gq(probes_betas(), title=title())
    print(plt)
})
output$cg_corrplot_gq = renderPlot(cg_corrplot_gq())

## Sample Analysis

plot_pca_sa = reactive({
    input$update_plot_pca_sa
    isolate({
        plt = .plot_pca_sa()
        print(plt)
    })})    
output$plot_pca_sa = renderPlot(plot_pca_sa())

conditionalPanel(condition="input.tabs == 'Probe Analysis'",
                 h5("Probe Analysis"),
                 wellPanel(
                     selectInput("analysis_select_gq", "select analysis",
                                 choices=c("", "probe summary", "probe correlation",
                                     "probe boxplot", "methylation track")),
                     uiOutput("select_cg_gq"),
                     conditionalPanel(condition=
                                      "input.analysis_select_gq == 'probe boxplot'",
                                      selectInput("comparison_select_gq", "comparison",
                                                  choices=c("tumor", "normal",
                                                      "tumor-normal",
                                                      "tumor-normal-stage"),
                                                  selected="tumor-normal")
                                      ),
                     uiOutput("cohort_select_gq"))
                 ),

conditionalPanel(condition="input.tabs == 'Sample Analysis'",
                 h5("Sample Analysis"),
                 wellPanel(
                     selectInput("analysis_select_sa", "select analysis",
                                 choices=c("PCA", "RPMM clusters")),

                     selectInput("probetypes_select_sa", "select probe location",
                                 choices=c("gene regulatory regions", "enhancers")),

                     selectInput("probetypes_select_sa", "select probe types",
                                 choices=c("CG Island", "CG shore")),
                     
                     sliderInput("cg_variance_sa", "probe variance percentile:",
                                 min=0, max=100, value=c(75,100), step=5),
                     uiOutput("cohort_select_sa"),
                     selectInput("samples_select_sa", "select samples",
                                 choices=c("all", "tumor", "normal")),
                     selectInput("comparison_select_sa", "select comparison",
                                 choices=c("tumor vs normal", "stage")),
                     
                     conditionalPanel(condition="input.analysis_select_sa == 'PCA'",
                                      actionButton("update_plot_pca_sa", "update PCA")),
                     
                     conditionalPanel(condition="input.analysis_select_sa == 'RPMM clusters'",
                                      actionButton("update_plot_rpmm_sa", "update RPMM"))
                     )
                 )

tabPanel("Probe Analysis",
         conditionalPanel(condition="input.analysis_select_gq == ''",
                          h5("select analysis, probes, and cohorts")),
         conditionalPanel(condition="input.analysis_select_gq == 'probe summary'",
                          h4("CpG info"),
                          tableOutput("cpg_info"),
                          br(),
                          h4("CpG islands info"),
                          tableOutput("cpg_islands_info")),                
         conditionalPanel(condition="input.analysis_select_gq == 'probe boxplot'",
                          h4("Methylation Boxplot"),
                          plotOutput("cg_boxplot_gq", width=800, height=600)),
         conditionalPanel(condition="input.analysis_select_gq == 'probe correlation'",
                          h4("Probe Correlation"),
                          plotOutput("cg_corrplot_gq", width=500)),
         conditionalPanel(condition="input.analysis_select_gq == 'methylation track'",
                          h4("Methylation Track"),
                          plotOutput("cg_track_gq", width=800, height=400))
         ),

tabPanel("Sample Analysis",
         conditionalPanel(condition="input.analysis_select_sa == ''",
                          h5("select analysis, samples")),

         conditionalPanel(condition="input.analysis_select_sa == 'PCA'",
                          plotOutput("plot_pca_sa", width=800)
                          ), 
         
         conditionalPanel(condition="input.analysis_select_sa == 'RPMM clusters'",
                          plotOutput("plot_rpmm_sa", width=800)
                          )

         ),

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
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressMessages(library("GenomicFeatures")) # GenomicRanges 


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
    ALIQUOT = aliquot[,c("bio.bcr_aliquot_uuid", "bio.bcr_sample_uuid", "admin.disease_code",
        "bio.sample_type_id")]
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

## Genome Annotation
## ENSEMBL_GTF = readRDS("tables/Homo_sapiens.GRCh37.73.rds")
## ENSEMBL_GTF = keepSeqlevels(ENSEMBL_GTF, c(1:22, "X", "Y"))
## ENSEMBL_GTF = renameSeqlevels(ENSEMBL_GTF, paste("chr", seqlevels(ENSEMBL_GTF), sep=""))

## HUGO2ENSG = data.table(as.data.frame(unique(mcols(ENSEMBL_GTF)[,c("gene_name", "gene_id")])))
## setkey(HUGO2ENSG, gene_name)
## ENSG2ENST = data.table(as.data.frame(unique(mcols(ENSEMBL_GTF)[,c("gene_id", "transcript_id")])))
## setkey(ENSG2ENST, gene_id)
