source(file.path(Sys.getenv("SYNC"), "libs", "lib.r"))

library(h5r)
library(data.table)

ifn = "tables/450k_tbl.h5"
ofn_row = "tables/450k_tbl_row.h5"

h5tbl = H5File(ifn, "r")
nodes = h5r::listH5Contents(h5tbl)
dsidx = laply(nodes, function(node) {
    node$type == 1
})
dspaths = names(nodes[dsidx])

h5_row = H5File(ofn_row, "w")
createH5Group(h5_row, "/450k")
createH5Group(h5_row, "/450k/TCGA")
for (dspath in dspaths) {
    print(dspath)
    ds = h5r::getH5Dataset(h5tbl, dspath)
    mx = ds[]
    ds_vars = as.matrix(rowVars(mx), ncol=1)
    createH5Group(h5_row, dspath)
    ds_vars_mx = createH5Dataset(h5_row, paste(dspath, "var", sep="/"),
        dims=dim(ds_vars), chunkSizes=c(dim(ds_vars)[[1]], 1), dType="double")
    writeH5Data(ds_vars_mx, ds_vars, c(1,1), c(dim(ds_vars)[[1]], dim(ds_vars)[[2]]))
    gc()
    gc()
    print("")
}



