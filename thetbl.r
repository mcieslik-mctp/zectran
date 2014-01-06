library(h5r)
library(stringr)

THETBL = H5File("tables/THETBL.h5", "r")
CONTENTS = readRDS("tables/THETBL.rds", "r")

.getChildren = function(root, depth=1) {
    tryCatch({
        root_name = ifelse(root@name==".", "/", paste("/", root@name, "/", sep=""))
        leafs = paste("/", names(CONTENTS)[2:length(CONTENTS)], sep="")
        matches = grepl(paste("^", root_name, sep=""), leafs)
        lnames = leafs[matches]
        rdepth = length(str_split(root_name, "/")[[1]])
        ldepth = sapply(str_split(lnames, "/"), length) - (rdepth - 1)
        values = lnames[ldepth %in% depth]
        values = str_replace(values, root_name, "")
    }, error=function(e) NULL)
}
.getH5Group = function(...) {
    tryCatch(getH5Group(...),
             error=function(e) NULL)
}
.getH5Dataset = function(...) {
    tryCatch(getH5Dataset(...),
             error=function(e) NULL)
}

getAnalyses = function(tbl=THETBL) {
    root = .getH5Group(tbl, ".")
    .getChildren(root)
}

getValues = function(analysis, tbl=THETBL) {
    root = .getH5Group(tbl, analysis)
    .getChildren(root)
}

getStudies = function(analysis, value, tbl=THETBL) {
    root = .getH5Group(tbl, paste(analysis, value, sep="/"))
    .getChildren(root)
}

getCohorts = function(analysis, value, study, tbl=THETBL) {
    root = .getH5Group(tbl, paste(analysis, value, study, sep="/"))
    .getChildren(root)
}

getVal = function(analysis, value, study, cohort, tbl=THETBL) {
    .getH5Dataset(tbl, paste(analysis, value, study, cohort, "val", sep="/"), inMemory=FALSE)    
}

getRid = function(analysis, value, study, cohort, tbl=THETBL) {
    .getH5Dataset(tbl, paste(analysis, value, study, cohort, "rid", sep="/"), inMemory=FALSE)[]
}

getCid = function(analysis, value, study, cohort, tbl=THETBL) {
    .getH5Dataset(tbl, paste(analysis, value, study, cohort, "cid", sep="/"), inMemory=FALSE)[]
}

getData = function(analysis, value, study, cohort, rows, cols, tbl=THETBL) {
    val = getVal(analysis, value, study, cohort, tbl)
    rid = getRid(analysis, value, study, cohort, tbl)
    cid = getCid(analysis, value, study, cohort, tbl)
    ridx = which(rid %in% rows)
    cidx = which(cid %in% cols)
    val[ridx, cidx]
}

