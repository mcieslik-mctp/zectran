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

shinyServer(function(input, output) {
    
    ## Gene Query
    
    output$cohort_select_gq = renderUI({
            selectInput("cohort_select_gq", label="select cohort(s)", choices=get_cohorts(), selected="PRAD",
                        multiple=TRUE)
    })

    output$cohort_select_sa = renderUI({
            selectInput("cohort_select_sa", label="select cohort(s)", choices=get_cohorts(), selected="PRAD",
                        multiple=TRUE)
    })

    
    ensgs = reactive({
        hgnc = input$query_hgnc_gq
        if (str_length(hgnc) > 0) {
            ensgs = hgnc2ensg(hgnc)$ensg
            if (length(ensgs) > 0) {
                return(ensgs)
            }
        }
        return(c("invalid gene"))
    })
    
    output$select_ensg_gq = renderUI(
        selectInput("select_ensg_gq", label="select gene", choices=ensgs())
        )
    
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
    # cg info
    cpg_info = reactive({
        .cpg_info(input$select_cg_gq)
    })
    output$cpg_info = renderTable(cpg_info())

    cpg_islands_info = reactive({
        .cpg_islands_info(input$select_cg_gq)
    })
    output$cpg_islands_info = renderTable(cpg_islands_info())

    # w/ beta values
    probes_betas = reactive({
        get_beta_n(input$cohort_select_gq, input$select_cg_gq)
    })
    title = reactive({
        sprintf("%s (%s) [%s]", input$query_hgnc_gq, input$select_ensg_gq, ensts_grl_select_region_str())
    })

    #
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
    

    ## Dataset Summary
    output$cohort_select_ds = renderUI({
        selectInput("cohort_select_ds", label="", choices=get_cohorts(), selected="PRAD",
                    multiple=TRUE)
    })
    # description
    cohort_description = reactive({.cohort_description(input$cohort_select_ds)})
    output$cohort_description = renderText(cohort_description())
    # summary
    cohort_summary = reactive({.cohort_summary(input$cohort_select_ds)})
    output$cohort_summary = renderTable(cohort_summary())
    # samples
    cohort_samples = reactive({.cohort_samples(input$cohort_select_ds)})
    output$cohort_samples = renderDataTable(cohort_samples())

    
})
print("server")
