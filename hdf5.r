get_cohorts = function() {
    cl = listH5Contents(H5FILE)
    tmp = cl[lapply(cl, "[[", "type") == 1]
    cohorts = unlist(lapply(tmp, "[", "name"), use.names=FALSE)
    return(cohorts)
}

get_samples = function(cohort) {
    ds = getH5Dataset(H5FILE, paste("/450k/TCGA", cohort, sep="/"), inMemory=FALSE)
    array_ids = getH5Attribute(ds, "ids")[]
    return(array_ids)
}

get_var = function(cohort, percentiles=FALSE) {
    ds = getH5Dataset(H5FILE_ROW, paste("/450k/TCGA", cohort, "var", sep="/"), inMemory=FALSE)
    mx = ds[]
    if (percentiles) {
        mx = as.integer(cut(mx, quantile(mx, probs=seq(0, 1, 0.01)), labels=1:100))
    }
    return(mx)
}

get_beta_1 = function(cohort, sample, cg) {
    cg_idx = match(cg, GR450K_IDS)
    ds = getH5Dataset(H5FILE, paste("/450k/TCGA", cohort, sep="/"), inMemory=FALSE)
    sm_idx = match(sample, getH5Attribute(ds, "ids")[])
    beta = ds[cg_idx, sm_idx]
    return(beta)
}

get_beta_n = function(cohorts, cgs) {
    if (length(cohorts) > 0 && (length(cgs) > 0)) {
        cgs_idx = which(GR450K_IDS %in% cgs)
        betas = llply(cohorts, function(cohort) {
            ds = getH5Dataset(H5FILE, paste("/450k/TCGA", cohort, sep="/"), inMemory=FALSE)
            ids = getH5Attribute(ds, "ids")[]
            cohort_betas = ds[cgs_idx,,]
            cohort_betas = matrix(cohort_betas, nrow=length(cgs), ncol=length(ids))
            colnames(cohort_betas) = ids
            rownames(cohort_betas) = cgs
            return(cohort_betas)
        })
        names(betas) = cohorts
        return(betas)
    }
}

print("hdf5")
