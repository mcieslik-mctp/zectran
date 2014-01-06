source("global.r")
source("thetbl.r")
source("tcga.r")
source("util.r")
source("convert.r")

H5DATA = "tables/450k_tbl.h5"
H5FILE = H5File(H5DATA, "r")
H5DATA_ROW = "tables/450k_tbl_row.h5"
H5FILE_ROW = H5File(H5DATA_ROW, "r")


