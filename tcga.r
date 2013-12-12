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

SAMPLE_COLS = c("bcr_sample_uuid","pathology_report_uuid","sample_type",
                "sample_type_id", "shortest_dimension")
getSamplesTable = function(samples, cohort, cols=SAMPLE_COLS) {
    cs = paste(cols, sep=",")
    query = sprintf(
        "select %s from biospecimen_sample_%s where lower(bcr_sample_uuid) in (?)", cs,
        tolower(cohort))
    unique(dbGetQuery(TCGA, query, samples))
}
getAllSamplesTable = function(cohort_samples, cols=SAMPLE_COLS) {
    tbls = list()
    for (cohort in names(cohort_samples)) {
        samples = cohort_samples[[cohort]]
        tbl = getSamplesTable(samples, cohort, cols)
        tbls[[cohort]] = tbl
    }
    rbind.fill(tbls)
}

getClinPatientJoin = function(patients, cohort, cols=PATIENT_COLS) {
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

filterAliquotsBySample = function(aliquots, cohort, column, select) {
    ch = tolower(cohort)
    map = as.data.frame(TCGAMAP[aliquots, bio.bcr_sample_uuid])
    q = dbGetQuery(TCGA, str_replace_all(sprintf(
    "select distinct bcr_sample_uuid,%s from biospecimen_sample_@ where bcr_sample_uuid in (?)",
        column), "@", ch), map$bio.bcr_sample_uuid)
    tbl = merge(x=map, y=q, by.x="bio.bcr_sample_uuid", by.y="bcr_sample_uuid")
    sel = (tbl[[column]] %in% select)
    sel_aliquots = as.character(tbl[sel,][["bio.bcr_aliquot_uuid"]])
}


filterAliquotsByPatient = function(aliquots, cohort, column, select) {
}






PATIENT_COLS = c("*")
getClinPatient = function(uuids, cohort="PRAD", cols=PATIENT_COLS) {
    cs = paste(cols, sep=",")
    query = sprintf(
        "select %s from clinical_patient_%s where lower(bcr_patient_uuid) in (?)", cs,
        tolower(cohort))
    dbGetQuery(TCGA, query, uuids)
}


cohortsToTable = function(analysis, value, study, cohorts, table_type) {
    tbls = list()
    for (cohort in cohorts) {
        aliquots = getCid(analysis, value, study, cohort)
        func = switch(table_type, "sample"={getBioSample}, "patient"={getBioPatient})
        col = switch(table_type, "sample"={"bio.bcr_sample_uuid"},
                                 "patient"={"shared.bcr_patient_uuid"})
        conv = as.character(TCGAMAP[aliquots, col][[col]])
        tbls[[cohort]] = func(samples, cohort)
    }
    tbl = rbind.fill(tbls)
    return(tbl)
}

## uniqueColVals = function(cohorts, table, column) {
##     vals = c()
##     for (cohort in cohorts) {
##         full_table = paste(table, tolower(cohort), sep="_")
##         query = sprintf("select distinct %s from %s", column, full_table)
##         vals = unique(c(vals, unique(dbGetQuery(TCGA, query)[[column]])))
##     }
##     return(vals)
## }
