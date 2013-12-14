library(stringr)
library(data.table)
library(RSQLite)

TCGAMAP = readRDS("tables/TCGAMAP.rds")
setkey(TCGAMAP, "bio.bcr_aliquot_uuid")

drv = dbDriver("SQLite")
TCGA = dbConnect(drv, "tables/clinical.sqlite")


convertAliquots = function(aliquots, to) {
    col = switch(to, "sample"={"bio.bcr_sample_uuid"},
                     "patient"={"shared.bcr_patient_uuid"})
    as.character(TCGAMAP[aliquots, get(col)]$V1)
}
convertAllAliquots = function(cohort_aliquots, to) {
    cohort_converts = list()
    for (cohort in names(cohort_aliquots)) {
        aliquots = cohort_aliquots[[cohort]]
        cohort_converts[[cohort]] = convertAliquots(aliquots, to)
    }
    return(cohort_converts)
}

getAliquots = function(analysis, value, study, cohort) {    
    aliquots = getCid(analysis, value, study, cohort)
}
getAllAliquots = function(analysis, value, study, cohorts) {
    cohort_aliquots = list()
    for (cohort in cohorts) {
        print(cohort)
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
    dbGetQuery(TCGA, query, samples)
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
filterAliquotsBySample = function(aliquots, cohort, column, select) {
    ch = tolower(cohort)
    map = as.data.frame(TCGAMAP[aliquots, bio.bcr_sample_uuid])
    q = dbGetQuery(TCGA, str_replace_all(sprintf(
    "select distinct bcr_sample_uuid,%s from biospecimen_sample_@ where lower(bcr_sample_uuid) in (?)",
        column), "@", ch), map$bio.bcr_sample_uuid)
    tbl = merge(x=map, y=q, by.x="bio.bcr_sample_uuid", by.y="bcr_sample_uuid")
    sel = (tbl[[column]] %in% select)
    sel_aliquots = as.character(tbl[sel,][["bio.bcr_aliquot_uuid"]])
}

filterAllAliquotsByFunction = function(cohort_aliquots, func, ...) {
    for (cohort in names(cohort_aliquots)) {
        aliquots = cohort_aliquots[[cohort]]
        aliquots_samplefiltered = filterAliquotsBySample(aliquots, cohort,
                input$samplecolumn_a_query_cb,
                input$samplefilter_a_select_cb,
                samplecolumn, samplefilters)

    }
}

##
getAllPatients = function(cohort) {
    query = sprintf("select distinct bcr_patient_uuid from clinical_patient_%s", tolower(cohort))
    dbGetQuery(TCGA, query)
}

PATIENT_COLS = c("*")
getPatientsTable = function(patients, cohort="PRAD", cols=PATIENT_COLS) {
    cs = paste(cols, sep=",")
    query = sprintf(
        "select distinct %s from clinical_patient_%s where lower(bcr_patient_uuid) in (?)", cs,
        tolower(cohort))
    dbGetQuery(TCGA, query, patients)
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
filterAliquotsByPatient = function(aliquots, cohort, column, select) {
    map = as.data.frame(TCGAMAP[aliquots, shared.bcr_patient_uuid])
    q = dbGetQuery(TCGA, str_replace_all(sprintf(
    "select distinct bcr_patient_uuid,%s from clinical_patient_@ where lower(bcr_patient_uuid) in (?)",
        column), "@", tolower(cohort)), map$shared.bcr_patient_uuid)
    tbl = merge(x=map, y=q, by.x="shared.bcr_patient_uuid", by.y="bcr_patient_uuid")
    sel = (tbl[[column]] %in% select)
    sel_aliquots = as.character(tbl[sel,][["bio.bcr_aliquot_uuid"]])
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
##         conv = as.character(TCGAMAP[aliquots, col][[col]])
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
