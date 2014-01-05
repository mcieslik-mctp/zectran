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
    
    ## ## Gene Query
    
    ## output$cohort_select_gq = renderUI({
    ##         selectInput("cohort_select_gq", label="select cohort(s)", choices=get_cohorts(), selected="PRAD",
    ##                     multiple=TRUE)
    ## })

    ## output$cohort_select_sa = renderUI({
    ##         selectInput("cohort_select_sa", label="select cohort(s)", choices=get_cohorts(), selected="PRAD",
    ##                     multiple=TRUE)
    ## })

    
    ## ensgs = reactive({
    ##     hgnc = input$query_hgnc_gq
    ##     if (str_length(hgnc) > 0) {
    ##         ensgs = hgnc2ensg(hgnc)$ensg
    ##         if (length(ensgs) > 0) {
    ##             return(ensgs)
    ##         }
    ##     }
    ##     return(c("invalid gene"))
    ## })
    
    ## output$select_ensg_gq = renderUI(
    ##     selectInput("select_ensg_gq", label="select gene", choices=ensgs())
    ##     )
    
    ## ensts = reactive({
    ##     if (is.null(input$select_ensg_gq) || input$select_ensg_gq == "invalid gene") {
    ##         return(NULL)
    ##     }
    ##     ensg = input$select_ensg_gq
    ##     ensts = ensg2enst(ensg)$enst
    ##     return(ensts)
    ## })
    
    ## output$select_enst_gq = renderUI(
    ##     selectInput("select_enst_gq", label="select transcript", choices=c("union", "longest", ensts()))
    ##     )

    ## ensts_grl = reactive({.ensts_grl(ensts())})

    ## ensts_tracks_gq = reactive({
    ##     input$plot_ensg_gq
    ##     isolate({
    ##         tks = .ensts_tracks_gq(ensts_grl())
    ##         print(tks)
    ##     })})
    ## output$ensts_tracks_gq = renderPlot(ensts_tracks_gq())
    
    
    ## ensts_grl_select = reactive({.ensts_grl_select(ensts_grl(), input$select_enst_gq)})

    ## ensts_grl_select_region = reactive({.ensts_grl_select_region(ensts_grl_select(), input$region_select_gq)})

    ## ensts_grl_select_region_str = reactive({.ensts_grl_select_region_str(ensts_grl_select_region())})

    ## output$query_region_limit_gq = renderUI(
    ##     textInput("query_region_limit_gq", label="limit region", ensts_grl_select_region_str()))

    ## ensts_grl_select_region_limit = reactive({.ensts_grl_select_region_str_parse(input$query_region_limit_gq)})

    ## ensts_grl_select_region_limit_probes = reactive({grl2probes(ensts_grl_select_region_limit())})

    ## ensts_select_tracks_gq = reactive({
    ##     input$plot_enst_gq
    ##     isolate({
    ##         tks = .ensts_select_tracks_gq(
    ##             ensts_grl_select(),
    ##             ensts_grl_select_region_limit(),
    ##             ensts_grl_select_region_limit_probes())
    ##         print(tks)
    ##     })})
    ## output$ensts_select_tracks_gq = renderPlot(ensts_select_tracks_gq())
    
    ## ensts_grl_select_region_limit_probes_names = reactive({.ensts_grl_select_region_limit_probes_names(
    ##     ensts_grl_select_region_limit_probes())})
    
    ## output$select_cg_gq = renderUI(selectInput("select_cg_gq", label="select probes",
    ##     choices=ensts_grl_select_region_limit_probes_names(), multiple=TRUE))

    ## ## Analyses
    ## # cg info
    ## cpg_info = reactive({
    ##     .cpg_info(input$select_cg_gq)
    ## })
    ## output$cpg_info = renderTable(cpg_info())

    ## cpg_islands_info = reactive({
    ##     .cpg_islands_info(input$select_cg_gq)
    ## })
    ## output$cpg_islands_info = renderTable(cpg_islands_info())

    ## # w/ beta values
    ## probes_betas = reactive({
    ##     get_beta_n(input$cohort_select_gq, input$select_cg_gq)
    ## })
    ## title = reactive({
    ##     sprintf("%s (%s) [%s]", input$query_hgnc_gq, input$select_ensg_gq, ensts_grl_select_region_str())
    ## })

    ## #
    ## cg_boxplot_gq = reactive({
    ##     plt = .cg_boxplot_gq(probes_betas(), comparison=input$comparison_select_gq, title=title())
    ##     print(plt)
    ## })
    ## output$cg_boxplot_gq = renderPlot(cg_boxplot_gq())
    
    ## cg_corrplot_gq = reactive({
    ##     plt = .cg_corrplot_gq(probes_betas(), title=title())
    ##     print(plt)
    ## })
    ## output$cg_corrplot_gq = renderPlot(cg_corrplot_gq())

    ## ## Sample Analysis

    ## plot_pca_sa = reactive({
    ##     input$update_plot_pca_sa
    ##     isolate({
    ##         plt = .plot_pca_sa()
    ##         print(plt)
    ##     })})    
    ## output$plot_pca_sa = renderPlot(plot_pca_sa())


    ## Dataset Summary
    all_cohorts = reactive({
        getCohorts(input$analysis_select_ds, input$value_select_ds, input$study_select_ds)
    })

    cohort_aliquots = reactive({
        cohort_aliquots = getAllAliquots(input$analysis_select_ds, input$value_select_ds, input$study_select_ds,
                       all_cohorts())
    })

    output$analysis_select_ds = renderUI({
        selectInput("analysis_select_ds", label="select analysis", choices=getAnalyses())
    })
    
    output$value_select_ds = renderUI({
        values = getValues(input$analysis_select_ds)
        if (is.null(values)) {
            return(NULL)
        }
        selectInput("value_select_ds", label="select value",
                    choices=values)
    })
    output$study_select_ds = renderUI({
        studies = getStudies(input$analysis_select_ds, input$value_select_ds)
        if (is.null(studies)) {
            return(NULL)
        }
        selectInput("study_select_ds", label="select study",
                    choices=studies)
    })
    output$dataset_table_ds = renderTable({
        if (is.null(cohort_aliquots())) {
            return(NULL)
        }
        ldply(names(cohort_aliquots()), function(cohort) {
            aliquots = cohort_aliquots()[[cohort]]
            n_aliquots = length(aliquots)
            n_tumor = length(unique(
                convertUUIDs(aliquots, "aliquot", "patient")))
            n_cohort_clinical = tryCatch({
                nrow(getAllPatients(cohort))
            }, error = function(e) {0L})
            df = data.frame(cohort, aliquots=length(aliquots), n_tumor, n_cohort_clinical)
            names(df) = c("cohort", "aliquots", "primary tumor", "TCGA cohort size (all samples w/ clinical data)")
            return(df)
        })
    })
    
    #### Cohort Builder
    cohort_samples = reactive({
        convertAllUUIDs(cohort_aliquots(), "aliquot", "sample")
    })
    
    cohort_patients = reactive({
        convertAllUUIDs(cohort_aliquots(), "aliquot", "patient")
    })

    ## CB - data - A
    a_cohort_aliquots = reactive({
        cohort_aliquots()[input$cohort_a_select_cb]
    })
    a_cohort_samples = reactive({
        cohort_samples()[input$cohort_a_select_cb]
    })
    a_cohort_patients = reactive({
        cohort_patients()[input$cohort_a_select_cb]
    })
    a_filtered_cohort_aliquots = reactive({
        filterAllAliquots(
            a_cohort_aliquots(),
            input$samplecolumn_a_query_cb,
            input$samplefilter_a_select_cb,
            input$patientcolumn_a_query_cb,
            input$patientfilter_a_select_cb)
    })
    a_filtered_cohort_patients = reactive({
        cohort_aliquots = a_filtered_cohort_aliquots()
        cohort_patients = lapply(convertAllUUIDs(cohort_aliquots, "aliquot", "patient"), unique)
    })
    a_refiltered_aliquots = reactive({
        convertUUIDs(input$patientselect_a_select_cb, "patient", "aliquot")
    })
    ## CB - ui - A
    output$cohort_a_select_cb = renderUI({
        cohorts = all_cohorts()
        if (is.null(cohorts)) {
            cohorts = "select dataset"
        }
        selectInput("cohort_a_select_cb", label="select cohort(s)", choices=cohorts,
                    multiple=TRUE)
    })    
    output$samplefilter_a_select_cb = renderUI({
        vals = getAllSamplesTable(a_cohort_samples(),
            input$samplecolumn_a_query_cb)[[input$samplecolumn_a_query_cb]]
        if (is.null(vals)) {
            vals = "select valid cohort(s) first"
        }
        lab = paste("select", input$samplecolumn_a_query_cb)
        selectInput("samplefilter_a_select_cb", label=lab, choices=sort(vals), multiple=TRUE)
    })
    output$patientfilter_a_select_cb = renderUI({
        vals = getAllPatientsTable(a_cohort_patients(),
            input$patientcolumn_a_query_cb)[[input$patientcolumn_a_query_cb]]
        if (is.null(vals)) {
            vals = "select valid cohort(s) first"
        }
        lab = paste("select", input$patientcolumn_a_query_cb)
        selectInput("patientfilter_a_select_cb", label=lab, choices=sort(vals),
                    multiple=TRUE)
    })
    output$patientselect_a_select_cb = renderUI({
        patients = as.vector(c(a_filtered_cohort_patients(), recursive=TRUE))
        selected = patients 
        if (length(patients) == 0) {
            patients = "select filter(s) first"
            selected = NULL
        }
        selectInput("patientselect_a_select_cb", label="select patient(s)", choices=patients,
                    selected=selected, multiple=TRUE)
    })

    ## CB - data - A
    b_cohort_aliquots = reactive({
        cohort_aliquots()[input$cohort_b_select_cb]
    })
    b_cohort_samples = reactive({
        cohort_samples()[input$cohort_b_select_cb]
    })
    b_cohort_patients = reactive({
        cohort_patients()[input$cohort_b_select_cb]
    })
    b_filtered_cohort_aliquots = reactive({
        filterAllAliquots(
            b_cohort_aliquots(),
            input$samplecolumn_b_query_cb,
            input$samplefilter_b_select_cb,
            input$patientcolumn_b_query_cb,
            input$patientfilter_b_select_cb)
    })
    b_filtered_cohort_patients = reactive({
        cohort_aliquots = b_filtered_cohort_aliquots()
        cohort_patients = lapply(convertAllUUIDs(cohort_aliquots, "aliquot", "patient"), unique)
    })
    b_refiltered_aliquots = reactive({
        convertUUIDs(input$patientselect_b_select_cb, "patient", "aliquot")
    })
    ## CB - ui - A
    output$cohort_b_select_cb = renderUI({
        cohorts = all_cohorts()
        if (is.null(cohorts)) {
            cohorts = "select dataset"
        }
        selectInput("cohort_b_select_cb", label="select cohort(s)", choices=cohorts,
                    multiple=TRUE)
    })    
    output$samplefilter_b_select_cb = renderUI({
        vals = getAllSamplesTable(b_cohort_samples(),
            input$samplecolumn_b_query_cb)[[input$samplecolumn_b_query_cb]]
        if (is.null(vals)) {
            vals = "select valid cohort(s) first"
        }
        lab = paste("select", input$samplecolumn_b_query_cb)
        selectInput("samplefilter_b_select_cb", label=lab, choices=sort(vals), multiple=TRUE)
    })
    output$patientfilter_b_select_cb = renderUI({
        vals = getAllPatientsTable(b_cohort_patients(),
            input$patientcolumn_b_query_cb)[[input$patientcolumn_b_query_cb]]
        if (is.null(vals)) {
            vals = "select valid cohort(s) first"
        }
        lab = paste("select", input$patientcolumn_b_query_cb)
        selectInput("patientfilter_b_select_cb", label=lab, choices=sort(vals),
                    multiple=TRUE)
    })
    output$patientselect_b_select_cb = renderUI({
        patients = as.vector(c(b_filtered_cohort_patients(), recursive=TRUE))
        selected = patients 
        if (length(patients) == 0) {
            patients = "select filter(s) first"
            selected = NULL
        }
        selectInput("patientselect_b_select_cb", label="select patient(s)", choices=patients,
                    selected=selected, multiple=TRUE)
    })

    

    ## ## CB - main   
    ## output$cohort_table_cb = renderTable({
    ##     input$cohorts_ab_action_cb
    ##     print("X")
    ##     isolate({
    ##         if (is.null(a_refiltered_aliquots()) || is.null(b_refiltered_aliquots())) {
    ##             return(NULL)
    ##         }
    ##         print("x")
    ##         print(a_refiltered_aliquots())
    ##         a_aic = length(a_refiltered_aliquots())
    ##         b_aic = length(b_refiltered_aliquots())           
    ##         df = data.frame(
    ##             c(a_aic, b_aic),
    ##             row.names=c("cohort A", "cohort B"))
    ##         colnames(df) = c("number of aliquots")
    ##         return(df)
    ##     })
    ## })


    ## ## CB - data - AB
    ## ab_filtered_cohort_aliquots = reactive({
    ##     ab_filters = list(
    ##         a=list(
    ##             cohort_aliquots=a_cohort_aliquots(),
    ##             samplecolumn=input$samplecolumn_a_query_cb,
    ##             samplefilters=input$samplefilter_a_select_cb,
    ##             patientcolumn=input$patientcolumn_a_query_cb,
    ##             patientfilters=input$patientfilter_a_select_cb),
    ##         b=list(
    ##             cohort_aliquots=b_cohort_aliquots(),
    ##             samplecolumn=input$samplecolumn_b_query_cb,
    ##             samplefilters=input$samplefilter_b_select_cb,
    ##             patientcolumn=input$patientcolumn_b_query_cb,
    ##             patientfilters=input$patientfilter_b_select_cb)
    ##         )
    ##     res=list()
    ##     for (ab in names(ab_filters)) {
    ##         res[[ab]] = list()            
    ##         f = ab_filters[[ab]]
    ##         cohort_aliquots = f[["cohort_aliquots"]]
    ##         samplecolumn = f[["samplecolumn"]]
    ##         samplefilters = f[["samplefilters"]]
    ##         patientcolumn = f[["patientcolumn"]]
    ##         patientfilters = f[["patientfilters"]]
    ##         for (cohort in names(cohort_aliquots)) {
    ##             aliquots = cohort_aliquots[[cohort]]
    ##             tryCatch({
    ##                 ## here aliquots are being filtered by patient and sample selects
    ##                 aliquots = filterAliquotsBySample(aliquots, cohort, samplecolumn,
    ##                     samplefilters)
    ##                 aliquots = filterAliquotsByPatient(aliquots, cohort, patientcolumn,
    ##                     patientfilters)
    ##             }, error=function(e) NULL)
    ##             res[[ab]][[cohort]] = aliquots
    ##         }
    ##     }
    ##     return(res)        
    ## })

    ## ab_filtered_cohort_patients = reactive({
    ##     inp = ab_filtered_cohort_aliquots()
    ##     res = list()
    ##     for (ab in names(inp)) {
    ##         cohort_aliquots = inp[[ab]]
    ##         ## multiple aliquots for the same patient need unique
    ##         cohort_patients = lapply(convertAllAliquots(cohort_aliquots, "patient"), unique)
    ##         res[[ab]] = cohort_patients
    ##     }
    ##     return(res)
    ## })
    ## patient_clinical_cb = reactive({
    ##     if (is.null(ab_filtered_cohort_patients())) {
    ##         return(NULL)
    ##     }
    ##     print("zz")
    ##     tbls = list()
    ##     for (ab in c("a", "b")) {
    ##         print(ab)
    ##         filtered_cohort_patients = ab_filtered_cohort_patients()[[ab]]
    ##         selected_patients = switch(ab,"a"={input$patient_a_select_cb},
    ##                                       "b"={input$patient_b_select_cb})
    ##         print(selected_patients)
    ##         for (cohort in names(filtered_cohort_patients)) {
    ##             print(cohort)
    ##             filtered_patients = filtered_cohort_patients[[cohort]]
    ##             selected_by_cohort = filtered_patients[filtered_patients %in% selected_patients]
    ##             if (length(selected_by_cohort) != 0) {
    ##                 print("zzz")
    ##                 tbls[[paste(ab, cohort, sep="_")]] = getPatientsTableJoin(selected_by_cohort,
    ##                         cohort)
    ##             }
    ##         }
    ##     }
    ##     unique(mcolsrbind.fillrbind.fill(tbls))
    ## })
    
    ## output$patient_clinical_cb = renderDataTable({
    ##     patient_clinical_cb()
    ## })

    ##
    ## output$cohorttally_table_da = renderTable({
    ##     if (is.null(ab_filtered_cohort_aliquots())) {
    ##         return(NULL)
    ##     }
    ##     rows = list()
    ##     for (ab in c("a", "b")) {
    ##         cohort_aliquots = ab_filtered_cohort_aliquots()[[ab]]
    ##         for (cohort in names(cohort_aliquots)) {
    ##             rows[[paste("cohort", toupper(ab), cohort)]] = list(
    ##                     cases=as.integer(length(cohort_aliquots[[cohort]])))
    ##         }
    ##     }
    ##     splat(rbind)(rows)
    ## })    
})
        
print("server")

source("thetbl.r")
source("tcga.r")
source("global_2.r")
