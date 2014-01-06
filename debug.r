source("global_2.r")
source("convert_2.r")
source("util.r")
source("thetbl.r")

ucsc = hgnc2ucsc("ZEB2")$ucsc[1]
grl = ucsc2grl(ucsc)
print(grl)
