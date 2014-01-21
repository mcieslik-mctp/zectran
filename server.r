source("thetbl.r")
source("tcga.r")
source("util.r")
source("convert.r")

shinyServer(function(input, output) {

    
    #### Individual Dataset Summary
    all_cohorts = reactive({
        getCohorts(input$analysis_select_ids, input$value_select_ids, input$study_select_ids)
    })

    cohort_aliquots = reactive({
        cohort_aliquots = getAllAliquots(input$analysis_select_ids, input$value_select_ids, input$study_select_ids,
                       all_cohorts())
    })

    output$analysis_select_ids = renderUI({
        selectInput("analysis_select_ids", label="select analysis", choices=getAnalyses())
    })
    
    output$value_select_ids = renderUI({
        values = getValues(input$analysis_select_ids)
        if (is.null(values)) {
            return(NULL)
        }
        selectInput("value_select_ids", label="select value",
                    choices=values)
    })
    output$study_select_ids = renderUI({
        studies = getStudies(input$analysis_select_ids, input$value_select_ids)
        if (is.null(studies)) {
            return(NULL)
        }
        selectInput("study_select_ids", label="select study",
                    choices=studies)
    })
    output$dataset_table_ids = renderTable({
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

    
    #### Cohort Builder (CB)
    cohort_samples = reactive({
        convertAllUUIDs(cohort_aliquots(), "aliquot", "sample")
    })
    cohort_patients = reactive({
        convertAllUUIDs(cohort_aliquots(), "aliquot", "patient")
    })

    ## CB - data - A
     a_cohort_aliquots = reactive({
        cohort_aliquots()[input$cohort_a_select_ids]
    })
    a_cohort_samples = reactive({
        cohort_samples()[input$cohort_a_select_ids]
    })
    a_cohort_patients = reactive({
        cohort_patients()[input$cohort_a_select_ids]
    })
    a_filtered_cohort_aliquots = reactive({
        filterAllAliquots(
            a_cohort_aliquots(),
            input$samplecolumn_a_query_ids,
            input$samplefilter_a_select_ids,
            input$patientcolumn_a_query_ids,
            input$patientfilter_a_select_ids)
    })
    a_filtered_cohort_patients = reactive({
        cohort_aliquots = a_filtered_cohort_aliquots()
        cohort_patients = lapply(convertAllUUIDs(cohort_aliquots, "aliquot", "patient"), unique)
    })
    a_refiltered_aliquots = reactive({
        a_filtered_aliquots = c(a_filtered_cohort_aliquots(), recursive=TRUE)
        a_patientrefiltered_aliquots = convertUUIDs(input$patientselect_a_select_ids, "patient", "aliquot")
        a_patientrefiltered_aliquots[a_patientrefiltered_aliquots %in% a_filtered_aliquots]
    })

    ## CB - ui - A
    output$cohort_a_select_ids = renderUI({
        cohorts = all_cohorts()
        if (is.null(cohorts)) {
            cohorts = "select dataset"
        }
        selectInput("cohort_a_select_ids", label="select cohort(s)", choices=cohorts,
                    multiple=TRUE)
    })    
    output$samplefilter_a_select_ids = renderUI({
        vals = getAllSamplesTable(a_cohort_samples(),
            input$samplecolumn_a_query_ids)[[input$samplecolumn_a_query_ids]]
        if (is.null(vals)) {
            vals = "select valid cohort(s) first"
        }
        lab = paste("select", input$samplecolumn_a_query_ids)
        selectInput("samplefilter_a_select_ids", label=lab, choices=sort(vals), multiple=TRUE)
    })
    output$patientfilter_a_select_ids = renderUI({
        vals = getAllPatientsTable(a_cohort_patients(),
            input$patientcolumn_a_query_ids)[[input$patientcolumn_a_query_ids]]
        if (is.null(vals)) {
            vals = "select valid cohort(s) first"
        }
        lab = paste("select", input$patientcolumn_a_query_ids)
        selectInput("patientfilter_a_select_ids", label=lab, choices=sort(vals),
                    multiple=TRUE)
    })
    output$patientselect_a_select_ids = renderUI({
        patients = as.vector(c(a_filtered_cohort_patients(), recursive=TRUE))
        selected = patients 
        if (length(patients) == 0) {
            patients = "select filter(s) first"
            selected = NULL
        }
        selectInput("patientselect_a_select_ids", label="select patient(s)", choices=patients,
                    selected=selected, multiple=TRUE)
    })

    ## CB - data - B
    b_cohort_aliquots = reactive({
        cohort_aliquots()[input$cohort_b_select_ids]
    })
    b_cohort_samples = reactive({
        cohort_samples()[input$cohort_b_select_ids]
    })
    b_cohort_patients = reactive({
        cohort_patients()[input$cohort_b_select_ids]
    })
    b_filtered_cohort_aliquots = reactive({
        filterAllAliquots(
            b_cohort_aliquots(),
            input$samplecolumn_b_query_ids,
            input$samplefilter_b_select_ids,
            input$patientcolumn_b_query_ids,
            input$patientfilter_b_select_ids)
    })
    b_filtered_cohort_patients = reactive({
        cohort_aliquots = b_filtered_cohort_aliquots()
        cohort_patients = lapply(convertAllUUIDs(cohort_aliquots, "aliquot", "patient"), unique)
    })
    b_refiltered_aliquots = reactive({
        b_filtered_aliquots = c(b_filtered_cohort_aliquots(), recursive=TRUE)
        b_patientrefiltered_aliquots = convertUUIDs(input$patientselect_b_select_ids, "patient", "aliquot")
        b_patientrefiltered_aliquots[b_patientrefiltered_aliquots %in% b_filtered_aliquots]
    })
    
    ## CB - ui - B
    output$cohort_b_select_ids = renderUI({
        cohorts = all_cohorts()
        if (is.null(cohorts)) {
            cohorts = "select dataset"
        }
        selectInput("cohort_b_select_ids", label="select cohort(s)", choices=cohorts,
                    multiple=TRUE)
    })    
    output$samplefilter_b_select_ids = renderUI({
        vals = getAllSamplesTable(b_cohort_samples(),
            input$samplecolumn_b_query_ids)[[input$samplecolumn_b_query_ids]]
        if (is.null(vals)) {
            vals = "select valid cohort(s) first"
        }
        lab = paste("select", input$samplecolumn_b_query_ids)
        selectInput("samplefilter_b_select_ids", label=lab, choices=sort(vals), multiple=TRUE)
    })
    output$patientfilter_b_select_ids = renderUI({
        vals = getAllPatientsTable(b_cohort_patients(),
            input$patientcolumn_b_query_ids)[[input$patientcolumn_b_query_ids]]
        if (is.null(vals)) {
            vals = "select valid cohort(s) first"
        }
        lab = paste("select", input$patientcolumn_b_query_ids)
        selectInput("patientfilter_b_select_ids", label=lab, choices=sort(vals),
                    multiple=TRUE)
    })
    output$patientselect_b_select_ids = renderUI({
        patients = as.vector(c(b_filtered_cohort_patients(), recursive=TRUE))
        selected = patients 
        if (length(patients) == 0) {
            patients = "select filter(s) first"
            selected = NULL
        }
        selectInput("patientselect_b_select_ids", label="select patient(s)", choices=patients,
                    selected=selected, multiple=TRUE)
    })

    ## CB - main
    patient_clinical_ids = reactive({
        tbls = list()
        for (ab in c("a", "b")) {
            filtered_cohort_patients = switch(ab, "a"=a_filtered_cohort_patients(),
                                                  "b"=b_filtered_cohort_patients()
                )
            selected_patients = switch(ab, "a"=input$patientselect_a_select_ids,
                                           "b"=input$patientselect_b_select_ids)
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
    output$patient_clinical_ids = renderDataTable({
        input$cohorts_ab_action_ids
        isolate({
            patient_clinical_ids()
        })
    })
    output$cohort_table_ids = renderTable({
        input$cohorts_ab_action_ids
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


    #### Paired Dataset Summary

    ## A
    output$analysis_a_select_pds = renderUI({
        selectInput("analysis_a_select_pds", label="select analysis", choices=getAnalyses())
    })
    output$value_a_select_pds = renderUI({
        values = getValues(input$analysis_a_select_pds)
        if (is.null(values)) {
            return(NULL)
        }
        selectInput("value_a_select_pds", label="select value",
                    choices=values)
    })
    output$study_a_select_pds = renderUI({
        studies = getStudies(input$analysis_a_select_pds, input$value_a_select_pds)
        if (is.null(studies)) {
            return(NULL)
        }
        selectInput("study_a_select_pds", label="select study",
                    choices=studies)
    })

    ## B
    output$analysis_b_select_pds = renderUI({
        selectInput("analysis_b_select_pds", label="select analysis", choices=getAnalyses())
    })
    output$value_b_select_pds = renderUI({
        values = getValues(input$analysis_b_select_pds)
        if (is.null(values)) {
            return(NULL)
        }
        selectInput("value_b_select_pds", label="select value",
                    choices=values)
    })
    output$study_b_select_pds = renderUI({
        studies = getStudies(input$analysis_b_select_pds, input$value_b_select_pds)
        if (is.null(studies)) {
            return(NULL)
        }
        selectInput("study_b_select_pds", label="select study",
                    choices=studies)
    })
    all_cohorts_pds = reactive({
        a_cohorts = getCohorts(input$analysis_a_select_pds, input$value_a_select_pds, input$study_a_select_pds)
        b_cohorts = getCohorts(input$analysis_b_select_pds, input$value_b_select_pds, input$study_b_select_pds)
        intersect(a_cohorts, b_cohorts)
    })
    output$cohort_select_pds = renderUI({
        cohorts = all_cohorts_pds()
        if (is.null(cohorts)) {
            cohorts = "select paired datasets"
        }
        selectInput("cohort_select_pds", label="select cohort(s)", choices=cohorts,
                    multiple=TRUE)
    })

    
    #### Methylation Analysis (MA)

    ## Gene Query (MA)
    output$name_text_ma = renderText({
        HGNC = toupper(input$hgnc_query_ma)
        name = hgnc2name(HGNC)$name
        if (length(name) == 0) {
            name = sprintf("gene symbol '%s' not found", HGNC)
        }
        return(name)
    })
    ucscs_grl_ma = reactive({
        HGNC = toupper(input$hgnc_query_ma)
        ucsc = ucsc2ucsc(hgnc2ucsc(HGNC)$ucsc)
        grl = ucsc2grl(ucsc)
        return(grl)
        
    })
    output$ucscmodel_plot_ma = renderPlot({
        input$ucscmodel_action_ma
        isolate({
            plt = tryCatch({
                grl = ucscs_grl_ma()
                plt = silentExec(
                    tracks(
                        autoplot(grl, group.selfish=TRUE, gap.geom="arrow") + theme_bw() + scale_x_sequnit("kb")
                        )
                    )
            }, error=function(e) {NULL})
            if (length(plt) != 0) {
                print(plt)
            }
        })
    })
    output$ucsc_select_ma = renderUI({
        txs = tryCatch(c("canonical", "union", names(ucscs_grl_ma())),
            error=function(e) "no transcripts found")
        selectInput("ucsc_select_ma", label="select transcript", choices=txs)})

    ucscs_grl_select_ma = reactive({
        transcripts = ucscs_grl_ma()
        if (length(input$ucsc_select_ma) != 0) {
            if (input$ucsc_select_ma == "canonical") {
                tx = getCanonicalTranscript(transcripts)
            } else if (input$ucsc_select_ma == "canonical") {
                tx = getUnionTranscript(transcripts)
            } else {
                tx = transcripts[input$ucsc_select_ma]
            }
        } else {
            tx = getEmptyTranscript()
        }
        return(tx)
    })
    ucscs_grl_select_region_ma = reactive({
        transcript = getFirstTranscript(ucscs_grl_select_ma())
        exon_number = mcols(transcript)$exon_number
        x = switch(input$region_select_ma,
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
               })
        return(x)
    })
    output$refinedregion_query_ma = renderUI({
        txt = tryCatch({
            gr = ucscs_grl_select_region_ma()
            sprintf("%s:%s-%s", seqnames(gr), start(gr), end(gr))
        }, error = function(e) "select transcript & region")
        textInput("refinedregion_query_ma", label="limit region", txt)
    })
    
    ucscs_grl_select_region_limit_ma = reactive({
        stringToRange(input$refinedregion_query_ma)
    })
    ucscs_grl_select_region_limit_probes_ma = reactive({
        grl2probes(ucscs_grl_select_region_limit_ma())
    })
    output$cpg_select_ma = renderUI({
        probe_names = tryCatch(
                names(ucscs_grl_select_region_limit_probes_ma()),
            error=function(e) {
                "select transcript & region"
            })
        selectInput("cpg_select_ma", label="select probes(s)", choices=probe_names, selected=probe_names, multiple=TRUE)
        })
    ucscs_grl_select_region_limit_probes_select_ma = reactive({
        ucscs_grl_select_region_limit_probes_ma()[input$cpg_select_ma]
    })
    output$probes_plot_ma = renderPlot({
        input$probes_action_ma
        isolate({
            plt = tryCatch({
                transcript = ucscs_grl_select_ma()
                region = ucscs_grl_select_region_limit_ma()
                probes = ucscs_grl_select_region_limit_probes_select_ma()
                silentExec(
                    tracks(
                        autoplot(transcript, group.selfish=TRUE, gap.geom="arrow") + theme_bw() + scale_x_sequnit("kb"),
                        autoplot(region) + theme_bw() + scale_x_sequnit("kb"),
                        autoplot(probes) + theme_bw() + scale_x_sequnit("kb")
                        )    
                    )
            }, error = function(e) {
                ## print(e)
                ## print(region)
                ## print(probes)
                NULL
            })
            if (length(plt) != 0) {
                print(plt)
            }
        })
    })
     output$dataset_text_ma = renderText({
         HTML(paste(
             paste("analysis:", input$analysis_select_ids),
             paste("value:", input$value_select_ids),
             paste("study:", input$study_select_ids),
             sep="<br>"))
         })

    #### Expression Analysis

    ## 
    output$name_text_ea = renderText({
        HGNC = toupper(input$hgnc_query_ea)
        name = hgnc2name(HGNC)$name
        if (length(name) == 0) {
            name = sprintf("gene symbol '%s' not found", HGNC)
        }
        return(name)
    })
    ucscs_grl_ea = reactive({
        HGNC = toupper(input$hgnc_query_ea)
        ucsc = ucsc2ucsc(hgnc2ucsc(HGNC)$ucsc)
        grl = ucsc2grl(ucsc)
        return(grl)
    })
    output$ucsc_select_ea = renderUI({
        txs = tryCatch(c("canonical", names(ucscs_grl_ea())),
            error=function(e) "no transcripts found")
        selectInput("ucsc_select_ea", label="select transcript", choices=txs)})
    ucscs_grl_select_ea = reactive({
        transcripts = ucscs_grl_ea()
        if (length(input$ucsc_select_ea) != 0) {
            if (input$ucsc_select_ea == "canonical") {
                tx = getCanonicalTranscript(transcripts)
            } else {
                tx = transcripts[input$ucsc_select_ea]
            }
        } else {
            tx = getEmptyTranscript()
        }
        return(tx)
    })

    output$ucscmodel_plot_ea = renderPlot({
        if (input$ucscmodel_action_ea > 0) {
            isolate({
                plt = tryCatch({
                txs = ucscs_grl_ea()
                tx = ucscs_grl_select_ea()
                mcols(txs)$selected = names(txs) %in% names(tx)
                plt = silentExec(
                    tracks(
                        autoplot(txs, group.selfish=TRUE, gap.geom="arrow", aes(fill=selected, color=selected)) + theme_bw() +
                        scale_x_sequnit("kb") +
                        scale_fill_manual(values=c("black", "red")) +
                        scale_color_manual(values=c("black", "red"))
                        )
                    )
            }, error=function(e) {print(e); NULL})
                if (length(plt) != 0) {
                    print(plt)
                }
            })
        }
    })

     output$dataset_text_ea = renderText({
         HTML(paste(
             paste("analysis:", input$analysis_select_ids),
             paste("value:", input$value_select_ids),
             paste("study:", input$study_select_ids),
             sep="<br>"))
     })


     output$cohorts_text_ea = renderText({
         HTML(paste(
             paste("cohort A"),
             paste("&nbsp&nbsp&nbsp&nbspsize:", length(a_refiltered_aliquots())),
             paste("cohort B"),
             paste("&nbsp&nbsp&nbsp&nbspsize:", length(b_refiltered_aliquots())),
             sep="<br>"))
     })
})

print("server")
