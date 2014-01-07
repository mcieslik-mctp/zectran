## shinyUI(pageWithSidebar(

##     headerPanel("Zectran the explorer"),

##     sidebarPanel(

##         conditionalPanel(condition="input.tabs == 'Dataset Summary'",
##                          h5("Select Dataset(s)"),
##                          wellPanel(
##                              uiOutput("analysis_select_ds"),
##                              uiOutput("value_select_ds"),
##                              uiOutput("study_select_ds")
##                              )
##                          ),
        
##         conditionalPanel(condition="input.tabs == 'Cohort Builder'",
##                      tags$div(class = "row-fluid",
##                               tags$div(class = "span6",
##                                        h5("(sub)cohort A"),
##                                        wellPanel(
##                                            uiOutput("cohort_a_select_cb"),
##                                            textInput("samplecolumn_a_query_cb", "sample filter column",
##                                                      value="sample_type_id"),
##                                            uiOutput("samplefilter_a_select_cb"),
##                                            textInput("patientcolumn_a_query_cb", "patient filter column",
##                                                      value="gleason_score"),
##                                            uiOutput("patientfilter_a_select_cb"),
##                                            uiOutput("patientselect_a_select_cb")
##                                            )
##                                        ),
##                               tags$div(class = "span6",
##                                        h5("(sub)cohort B"),
##                                        wellPanel(
##                                            uiOutput("cohort_b_select_cb"),
##                                            textInput("samplecolumn_b_query_cb", "Sample Column",
##                                                      value="sample_type_id"),
##                                            uiOutput("samplefilter_b_select_cb"),
##                                            textInput("patientcolumn_b_query_cb", "patient filter column",
##                                                      value="gleason_score"),
##                                            uiOutput("patientfilter_b_select_cb"),
##                                            uiOutput("patientselect_b_select_cb")
##                                            )
##                                        )
##                               ),
##                          actionButton("cohorts_ab_action_cb", "update cohorts")
##                          ),
        
##         conditionalPanel(condition="input.tabs == 'Methylation Analysis'",
##                          h5("Select Gene"),
##                          wellPanel(
##                              tags$div(class = "row-fluid",
##                                       tags$div(class = "span6", textInput("hgnc_query_me", "HUGO (HGNC) gene name")),
##                                       tags$div(class = "span6", textOutput("name_text_me"))
##                              ),
##                              actionButton("ucscmodel_action_me", "plot gene model")
##                              ),
##                          h5("Refine Region"),
##                          wellPanel(
##                              uiOutput("ucsc_select_me"),
##                              selectInput("region_select_me", "select region",
##                                          choices=c(
##                                              "TSS slim (500, 200)",
##                                              "TSS wide (2000, 500)",
##                                              "gene regulatory region")),
##                              uiOutput("refinedregion_query_me"),
##                              actionButton("ucscregion_plot_me", "plot refined region")
##                              )
##                          )
##     ),
    
##     mainPanel(
        
##         tabsetPanel(

##             tabPanel("Dataset Summary",
##                      h4("Dataset Summary"),
##                      tableOutput("dataset_table_ds")
##                      ),

##             tabPanel("Cohort Builder",
##                      h4("Patient Inspector"),
##                      dataTableOutput("patient_clinical_cb"),
##                      h4("Cohort Summary"),
##                      tableOutput("cohort_table_cb")
##                      ),

##             tabPanel("Methylation Analysis",
##                      h4("UCSC Gene Model"),
##                      plotOutput("ucscmodel_plot_me", width=800),
##                      h4("Selected UCSC Transcript"),
##                      plotOutput("ucsctranscript_plot_me", width=800, height=300)
##                      ),
            
##             id="tabs")
##     )
## ))

## Tables






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

get_cohorts = function() {
    cl = listH5Contents(H5FILE)
    tmp = cl[lapply(cl, "[[", "type") == 1]
    cohorts = unlist(lapply(tmp, "[", "name"), use.names=FALSE)
    return(cohorts)
}

get_samples = function(cohort) {
    ds = getH5Dataset(H5FILE, paste("/450k/TCGA", cohort, sep="/"), inMemory=FALSE)
    array_ids = getH5Attribute(ds, "ids")[]
    return(array_ids)
}

get_var = function(cohort, percentiles=FALSE) {
    ds = getH5Dataset(H5FILE_ROW, paste("/450k/TCGA", cohort, "var", sep="/"), inMemory=FALSE)
    mx = ds[]
    if (percentiles) {
        mx = as.integer(cut(mx, quantile(mx, probs=seq(0, 1, 0.01)), labels=1:100))
    }
    return(mx)
}

get_beta_1 = function(cohort, sample, cg) {
    cg_idx = match(cg, GR450K_IDS)
    ds = getH5Dataset(H5FILE, paste("/450k/TCGA", cohort, sep="/"), inMemory=FALSE)
    sm_idx = match(sample, getH5Attribute(ds, "ids")[])
    beta = ds[cg_idx, sm_idx]
    return(beta)
}

get_beta_n = function(cohorts, cgs) {
    if (length(cohorts) > 0 && (length(cgs) > 0)) {
        cgs_idx = which(GR450K_IDS %in% cgs)
        betas = llply(cohorts, function(cohort) {
            ds = getH5Dataset(H5FILE, paste("/450k/TCGA", cohort, sep="/"), inMemory=FALSE)
            ids = getH5Attribute(ds, "ids")[]
            cohort_betas = ds[cgs_idx,,]
            cohort_betas = matrix(cohort_betas, nrow=length(cgs), ncol=length(ids))
            colnames(cohort_betas) = ids
            rownames(cohort_betas) = cgs
            return(cohort_betas)
        })
        names(betas) = cohorts
        return(betas)
    }
}

