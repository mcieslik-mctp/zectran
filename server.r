source("thetbl.r")
source("tcga.r")
source("util.r")
source("convert.r")

shinyServer(function(input, output) {
    
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
        a_filtered_aliquots = c(a_filtered_cohort_aliquots(), recursive=TRUE)
        a_patientrefiltered_aliquots = convertUUIDs(input$patientselect_a_select_cb, "patient", "aliquot")
        a_patientrefiltered_aliquots[a_patientrefiltered_aliquots %in% a_filtered_aliquots]
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

    ## CB - data - B
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
        b_filtered_aliquots = c(b_filtered_cohort_aliquots(), recursive=TRUE)
        b_patientrefiltered_aliquots = convertUUIDs(input$patientselect_b_select_cb, "patient", "aliquot")
        b_patientrefiltered_aliquots[b_patientrefiltered_aliquots %in% b_filtered_aliquots]
    })
    
    ## CB - ui - B
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

    ## CB - main
    patient_clinical_cb = reactive({
        tbls = list()
        for (ab in c("a", "b")) {
            filtered_cohort_patients = switch(ab, "a"=a_filtered_cohort_patients(),
                                                  "b"=b_filtered_cohort_patients()
                )
            selected_patients = switch(ab, "a"=input$patientselect_a_select_cb,
                                           "b"=input$patientselect_b_select_cb)
            for (cohort in names(filtered_cohort_patients)) {
                filtered_patients = filtered_cohort_patients[[cohort]]
                selected_by_cohort = filtered_patients[filtered_patients %in% selected_patients]
                if (length(selected_by_cohort) != 0) {
                    tbls[[paste(ab, cohort, sep="_")]] = getPatientsTableJoin(selected_by_cohort,
                            cohort)
                }
            }
        }
        unique(rbind.fill(tbls))
    })
    
    output$patient_clinical_cb = renderDataTable({
        input$cohorts_ab_action_cb
        isolate({
            patient_clinical_cb()
        })
    })
    
    output$cohort_table_cb = renderTable({
        input$cohorts_ab_action_cb
        isolate({
            a_aic = length(a_refiltered_aliquots())
            b_aic = length(b_refiltered_aliquots())            
            df = data.frame(
                c(a_aic, b_aic),
                row.names=c("cohort A", "cohort B"))
            colnames(df) = c("number of aliquots")
            return(df)
        })
    })

    output$name_text_me = renderText({
        HGNC = toupper(input$hgnc_query_me)
        name = hgnc2name(HGNC)$name
        if (length(name) == 0) {
            name = sprintf("gene symbol '%s' not found", HGNC)
        }
        return(name)
    })

    ucscs = reactive({
        HGNC = toupper(input$hgnc_query_me)
        hgnc2ucsc(HGNC)$ucsc        
    })

    ucscs_grl = reactive({
        ucsc2grl(ucscs())
    })
    
    ucscmodel_plot_me = reactive({
        input$ucscmodel_action_me
        isolate({
            grl = ucscs_grl()
            plt = tracks(
                autoplot(grl, group.selfish=TRUE, gap.geom="arrow") + theme_bw() + scale_x_sequnit("kb")
                )
            print(plt)
        })})
    
    output$ucscmodel_plot_me = renderPlot(ucscmodel_plot_me())

    ###
    
    ucsc_select_me = reactive({
        txs = ucscs()
        if (length(txs) == 0) {
            return("no transcripts found") 
        }
        txs = c("union", "longest", txs)
        return(txs)
    })
    
    output$ucsc_select_me = renderUI(
        selectInput("ucsc_select_me", label="select transcript", choices=ucsc_select_me()))
    
    ucscs_grl_select = reactive({
        grl = ucscs_grl()
        select = input$ucsc_select_me
        switch(select,
               "union"={
                   all = unlist(grl)
                   uni = reduce(all)
                   mcols(uni)$type = "exon"
                   mcols(uni)$exon_number = ifelse(strand(uni) == "+", 1:length(uni), length(uni):1)
                   unis = GRangesList()
                   gene_id = unique(mcols(all)$gene_id)
                   mcols(uni)$gene_id = gene_id
                   unis[[gene_id]] = uni
                   unis
               },
               "longest"={
                   grl_width = unlist(lapply(grl, function(transcript) {
                       max(end(transcript)) - min(start(transcript))
                   }))
                   select = names(grl_width[order(grl_width, decreasing=TRUE)[1]])
                   grl[select]
               },
               grl[select]
               )
    })

    ucscs_grl_select_region = reactive({
        transcript = ucscs_grl_select()[[1]] # should be only one
        exon_number = mcols(transcript)$exon_number
        switch(input$region_select_me,
               "TSS slim (500, 200)"={
                   first_exon = transcript[exon_number == 1]
                   promoters(first_exon, 500, 200)
               },
               "TSS wide (2000, 500)"={
                   first_exon = transcript[exon_number == 1]
                   promoters(first_exon, 2000, 500)
               },
               "gene regulatory region"={
                   flankTranscript(transcript, 2000, 2000)
               }
               )
    })
    
    ucscs_grl_select_region_str = reactive({
        print(ucscs_grl_select_region())
        sprintf("%s:%s-%s", seqnames(ucscs_grl_select_region()),
                start(ucscs_grl_select_region()),
                end(ucscs_grl_select_region()))
    })

    output$refinedregion_query_me = renderUI(
        textInput("refinedregion_query_me", label="limit region",
                  tryCatch(ucscs_grl_select_region_str(), error=function(e) "select transcript & region")))

    ucscs_grl_select_region_limit = reactive({
        input$refinedregion_query_me
        
        str_replace_all(input$refinedregion_query_me, ",", "")
        fields = unlist(str_split(str_split(input$refinedregion_query_me, pattern=":")[[1]], "-"))
        chr = fields[1]
        start = as.integer(fields[2])
        end = as.integer(fields[3])
        gr = GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), strand="*", type="range")
        return(gr)

    })
    

    



    







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



    

    
    

})
suppressMessages(library("ggbio"))

print("server")
