source("thetbl.r")
source("tcga.r")
source("util.r")
source("convert.r")

#source("convert_old.r")
#source("global_old.r")

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

    ####
    ucscs = reactive({
        hgnc = input$hgnc_query_gq
        if (str_length(hgnc) > 0) {
            ucscs = hgnc2ucsc(hgnc)$ucsc
            if (length(ucscs) > 0) {
                return(ucscs)
            }
        }
        return(c("invalid gene"))
    })
    
    output$ucsc_select_gq = renderUI(
        selectInput("select_ensg_gq", label="select gene", choices=ensgs()))
    

})

print("server")
