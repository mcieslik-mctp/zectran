library(stringr)
library(data.table)
library(RSQLite)

TCGAMAP_A = readRDS("tables/TCGAMAP.rds")
setkey(TCGAMAP_A, "bio.bcr_aliquot_uuid")
TCGAMAP_P =  copy(TCGAMAP_A)
setkey(TCGAMAP_P, "shared.bcr_patient_uuid")

drv = dbDriver("SQLite")
TCGA = dbConnect(drv, "tables/clinical.sqlite")

convertUUIDs = function(uuids, from, to) {
    tryCatch({
        MAP = switch(from, "aliquot"={TCGAMAP_A},
            "patient"={TCGAMAP_P})
        col = switch(to, "sample"={"bio.bcr_sample_uuid"},
            "aliquot"={"bio.bcr_aliquot_uuid"},
            "patient"={"shared.bcr_patient_uuid"})
        tbl = MAP[uuids, get(col)]
        as.character(tbl$V1)
    }, error=function(e) {
        #print(e)
        NULL
    })
}

convertAllUUIDs = function(cohort_uuids, from, to) {
    cohort_converts = list()
    for (cohort in names(cohort_uuids)) {
        uuids = cohort_uuids[[cohort]]
        cohort_converts[[cohort]] = convertUUIDs(uuids, from, to)
    }
    return(cohort_converts)
}

getAliquots = function(analysis, value, study, cohort) {    
    aliquots = getCid(analysis, value, study, cohort)
}
getAllAliquots = function(analysis, value, study, cohorts) {
    cohort_aliquots = list()
    for (cohort in cohorts) {
        cohort_aliquots[[cohort]] = getAliquots(analysis, value, study, cohort)
    }
    return(cohort_aliquots)
}

##
SAMPLE_COLS = c("bcr_sample_uuid","pathology_report_uuid","sample_type",
                "sample_type_id", "shortest_dimension")
getSamplesTable = function(samples, cohort, cols=SAMPLE_COLS) {
    cs = paste(cols, sep=",")
    query = sprintf(
        "select distinct %s from biospecimen_sample_%s where lower(bcr_sample_uuid) in (?)", cs,
        tolower(cohort))
    tryCatch(dbGetQuery(TCGA, query, samples), error=function(e) {
        #print(e)
        NULL
    })
}
getAllSamplesTable = function(cohort_samples, cols=SAMPLE_COLS) {
    tbls = list()
    for (cohort in names(cohort_samples)) {
        samples = cohort_samples[[cohort]]
        tbl = getSamplesTable(samples, cohort, cols)
        tbls[[cohort]] = tbl
    }
    unique(rbind.fill(tbls))
}
filterAliquotsBySampleFeatures = function(aliquots, cohort, column, select) {
    if (length(aliquots) == 0) {
        return(NULL)
    }
    sel_aliquots = tryCatch(
        {
            map = as.data.frame(TCGAMAP_A[aliquots, bio.bcr_sample_uuid])
            template = paste(
                "select distinct",
                "lower(bcr_sample_uuid) as bcr_sample_uuid, %s",
                "from biospecimen_sample_@ where lower(bcr_sample_uuid) in (?)")
            query = str_replace_all(sprintf(template, column), "@", tolower(cohort))
            res = dbGetQuery(TCGA, query, map$bio.bcr_sample_uuid)
            tbl = merge(x=map, y=res, by.x="bio.bcr_sample_uuid", by.y="bcr_sample_uuid")
            tbl = tbl[!duplicated(tbl),] # why oh why
            sel = (tbl[[column]] %in% select)
            sel_aliquots = as.character(tbl[sel,][["bio.bcr_aliquot_uuid"]])
        },
        error = function(e) {print(e); NULL})
    return(sel_aliquots)
}

##
getAllPatients = function(cohort) {
    query = sprintf("select distinct lower(bcr_patient_uuid) as bcr_patient_uuid from clinical_patient_%s",
        tolower(cohort))
    dbGetQuery(TCGA, query)
}

PATIENT_COLS = c("*")
getPatientsTable = function(patients, cohort="PRAD", cols=PATIENT_COLS) {
    cs = paste(cols, sep=",")
    query = sprintf(
        "select distinct %s from clinical_patient_%s where lower(bcr_patient_uuid) in (?)", cs,
        tolower(cohort))
    tryCatch(dbGetQuery(TCGA, query, patients), error=function(e) {
        #print(e)
        NULL
    })
}
getAllPatientsTable = function(cohort_patients, cols=PATIENT_COLS) {
    tbls = list()
    for (cohort in names(cohort_patients)) {
        patients = cohort_patients[[cohort]]
        tbl = getPatientsTable(patients, cohort, cols)
        tbls[[cohort]] = tbl
    }
    unique(rbind.fill(tbls))
}

filterAliquotsByPatientFeatures = function(aliquots, cohort, column, select) {
    if (length(aliquots) == 0) {
        return(NULL)
    }
    sel_aliquots = tryCatch(
        {
            map = as.data.frame(TCGAMAP_A[aliquots, shared.bcr_patient_uuid])
            template = paste(
                "select distinct",
                "lower(bcr_patient_uuid) as bcr_patient_uuid, %s",
                "from clinical_patient_@ where lower(bcr_patient_uuid) in (?)")
            query = str_replace_all(sprintf(template, column), "@", tolower(cohort))            
            res = dbGetQuery(TCGA, query, map$shared.bcr_patient_uuid)
            tbl = merge(x=map, y=res, by.x="shared.bcr_patient_uuid", by.y="bcr_patient_uuid")
            tbl = tbl[!duplicated(tbl),] # why oh why
            sel = (tbl[[column]] %in% select)
            as.character(tbl[sel,][["bio.bcr_aliquot_uuid"]])
        },
        error = function(e) {
            print(e)
            NULL}
        )
    return(sel_aliquots)
}


filterAllAliquots = function(cohort_aliquots, samplecolumn, samplefilters, patientcolumn, patientfilters) {
    res = list()
    for (cohort in names(cohort_aliquots)) {
        aliquots = cohort_aliquots[[cohort]]
        print("a")
        print(length(aliquots))
        aliquots = filterAliquotsBySampleFeatures(aliquots, cohort, samplecolumn,
            samplefilters)
        print(length(aliquots))
        aliquots = filterAliquotsByPatientFeatures(aliquots, cohort, patientcolumn,
            patientfilters)
        print(length(aliquots))
        print("z")
        res[[cohort]] = aliquots
    }
    return(res)
}

getPatientsTableJoin = function(patients, cohort, cols=PATIENT_COLS) {
    ch = tolower(cohort)
    tables = dbListTables(TCGA)
    has_patient = paste("clinical_patient", ch, sep="_") %in% tables
    if (!has_patient) {
        return(NULL)
    }
    has_cqcf = paste("clinical_cqcf", ch, sep="_") %in% tables
    has_drug = paste("clinical_drug", ch, sep="_") %in% tables
    has_radi = paste("clinical_radiation", ch, sep="_") %in% tables
    tbl = dbGetQuery(TCGA,
               str_replace_all(
                   paste(
                       str_replace("select distinct YYY from clinical_patient_XXX", "YYY",
                                   paste(cols,sep=",")),
                       ## cqcf
                       ifelse(has_cqcf,
                              paste(
                                  "left join clinical_cqcf_XXX on",
                                  "upper(clinical_cqcf_XXX.bcr_patient_barcode) =",
                                  "upper(clinical_patient_XXX.bcr_patient_barcode)"),
                              ""),
                       ## drug
                       ifelse(has_drug,
                              paste(
                                  "left join clinical_drug_XXX on",
                                  "upper(clinical_drug_XXX.bcr_patient_barcode) =",
                                  "upper(clinical_patient_XXX.bcr_patient_barcode)"),
                              ""),
                       ## radiation
                       ifelse(has_radi,
                              paste(
                                  "left join clinical_radiation_XXX on",
                                  "upper(clinical_radiation_XXX.bcr_patient_barcode) =",
                                  "upper(clinical_patient_XXX.bcr_patient_barcode)"),
                              ""),
                       ## where
                       "where lower(clinical_patient_XXX.bcr_patient_uuid) in (?)",
                       ""
                       ), "XXX", ch), patients)
    return(unique(tbl))
}




## cohortsToTable = function(analysis, value, study, cohorts, table_type) {
##     tbls = list()
##     for (cohort in cohorts) {
##         aliquots = getCid(analysis, value, study, cohort)
##         func = switch(table_type, "sample"={getBioSample}, "patient"={getBioPatient})
##         col = switch(table_type, "sample"={"bio.bcr_sample_uuid"},
##                                  "patient"={"shared.bcr_patient_uuid"})
##         convertAliquots = as.character(TCGAMAP[aliquots, col][[col]])
##         tbls[[cohort]] = func(samples, cohort)
##     }
##     tbl = rbind.fill(tbls)
##     return(tbl)
## }
## uniqueColVals = function(cohorts, table, column) {
##     vals = c()
##     for (cohort in cohorts) {
##         full_table = paste(table, tolower(cohort), sep="_")
##         query = sprintf("select distinct %s from %s", column, full_table)
##         vals = unique(c(vals, unique(dbGetQuery(TCGA, query)[[column]])))
##     }
##     return(vals)
## }

    ## # description
    ## cohort_description = reactive({.cohort_description(input$cohort_select_ds)})
    
    ## # summary
    ## cohort_summary = reactive({.cohort_summary(input$cohort_select_ds)})
    ## output$cohort_summary = renderTable(cohort_summary())
    ## # samples
    ## cohort_samples = reactive({.cohort_samples(input$cohort_select_ds)})
    ## #output$cohort_samples = renderDataTable(cohort_samples())
