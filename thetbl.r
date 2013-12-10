library(h5r)
library(stringr)

THETBL = H5File("tables/thetbl.h5", "r")

.getChildren = function(root) {
    leafs = listH5Contents(root)
    lnames = names(leafs)[2:length(leafs)]
    ldepth = sapply(str_split(lnames, "/"), length)
    values = lnames[ldepth == 1]
}

getAnalyses = function(tbl=THETBL) {
    .getChildren(tbl)
}

getValues = function(analysis, tbl=THETBL) {
    root = getH5Group(tbl, analysis)
    .getChildren(root)
}

getStudies = function(analysis, value, tbl=THETBL) {
    root = getH5Group(tbl, paste(analysis, value, sep="/"))
    .getChildren(root)
}

getCohorts = function(analysis, value, study, tbl=THETBL) {
    root = getH5Group(tbl, paste(analysis, value, study, sep="/"))
    .getChildren(root)
}

getTable = function(analysis, value, study, cohort, tbl=THETBL) {
    ds = getH5Dataset(tbl, paste(analysis, value, study, cohort, sep="/"), inMemory=FALSE)    
}

getSamples = function(analysis, value, study, cohort, tbl=THETBL) {
    ds = getH5Dataset(tbl, paste(analysis, value, study, cohort, sep="/"), inMemory=FALSE)
    samples = getH5Attribute(ds, attrName="ids")[]
}

getIndex = function(analysis, value, tbl=THETBL) {
    index = getH5Dataset(tbl, paste(analysis, value, "idx", sep="/"))[]
}

getData = function(analysis, value, study, cohort, samples, ids, tbl=THETBL) {
    table = getTable(analysis, value, study, cohort, tbl)
    all_samples = getSamples(analysis, value, study, cohort, tbl)
    index = getIndex(analysis, value, tbl)
    rows = which(index == ids)
    cols = which(all_samples == samples)
    table[rows, cols]
}
