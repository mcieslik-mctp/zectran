source("global.r")
source("thetbl.r")
source("tcga.r")
source("util.r")
source("convert.r")

tread = function(fn, sep="\t") {
    x = read.csv(fn, header=TRUE, sep=sep, comment.char="")
    first = str_replace(names(x)[1], "X.", "")
    tmp = if (first=="") "id" else first
    names(x) = c(tmp, names(x)[2:length(x)])
    attr(x, "first") = first
    return(x)
}

ALIQUOT_FN = "tables/bio:bcr_aliquot_uuid.tsv"
SDRF_FN = "tables/mage-450k_merged_sdrf.txt"
make_aliquot_meta = function(aliquot_fn=ALIQUOT_FN, sdrf_fn=SDRF_FN) {
    aliquot = tread(aliquot_fn)
    aliquot = aliquot[!is.na(aliquot$bio.bcr_aliquot_uuid),]
    ALIQUOT = aliquot[,c("bio.bcr_aliquot_uuid", "bio.bcr_sample_uuid", "admin.disease_code",
        "bio.sample_type_id")]
    ALIQUOT$bio.bcr_aliquot_uuid = toupper(ALIQUOT$bio.bcr_aliquot_uuid)
    ALIQUOT$bio.bcr_sample_uuid = toupper(ALIQUOT$bio.bcr_sample_uuid)
    
    sdrf = tread(sdrf_fn)
    ARRAY = sdrf[c("Extract.Name", "Array.Data.File")]
    names(ARRAY) = c("bio.bcr_aliquot_uuid", "file")
    ARRAY$bio.bcr_aliquot_uuid = toupper(ARRAY$bio.bcr_aliquot_uuid)

    ARRAY$array_450k_id = str_sub(ARRAY$file, end=-10)
    ARRAY$file = NULL
    ARRAY = unique(ARRAY)
    ALIMETA = merge(ALIQUOT, ARRAY, by="bio.bcr_aliquot_uuid")
}
ALIMETA = make_aliquot_meta()

## PRAD
cohort = "PRAD"
cgs = GR_450K_IDS[c(10,20,100,1000,10000,777,8321,31928)]

get_samples = function(cohort) {
    ds = getH5Dataset(H5FILE, paste("/450k/TCGA", cohort, sep="/"), inMemory=FALSE)
    array_ids = getH5Attribute(ds, "ids")[]
    return(array_ids)
}

sample = get_samples(cohort)[[20]]
aliquot_uuid = tolower(ALIMETA[ALIMETA$array_450k_id == sample, "bio.bcr_aliquot_uuid"])

## OLD
cgs_idx = which(GR450K_IDS %in% cgs)
H5DATA = "tables/450k_tbl.h5"
H5FILE = H5File(H5DATA, "r")
H5DATA_ROW = "tables/450k_tbl_row.h5"
H5FILE_ROW = H5File(H5DATA_ROW, "r")

get_beta_1 = function(cohort, sample, cg) {
    cg_idx = match(cg, GR_450K_IDS)
    ds = getH5Dataset(H5FILE, paste("/450k/TCGA", cohort, sep="/"), inMemory=FALSE)
    sm_idx = match(sample, getH5Attribute(ds, "ids")[])
    beta = ds[cg_idx, sm_idx]
    return(beta)
}

for (cg in cgs) {
    print("#")
    print(get_beta_1(cohort, sample, cg))
    print(getData("HumanMethylation450.Level_2", "lvl-2.beta", "TCGA", "PRAD", cg, aliquot_uuid))
}
