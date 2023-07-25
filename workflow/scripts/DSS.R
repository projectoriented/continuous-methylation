#!/usr/bin/env Rscript
# Author: Mei Wu, https://github.com/projectoriented

####### -------------- Libraries -------------- #######
library("DSS")

####### -------------- Functions -------------- #######
process_files <- function(file_names, sample_names) {
    # Create an empty list to store the data frames
    data_list <- list()

    # Iterate over the file names
    for (i in seq_along(file_names)) {
        # Read the data from each file
        data <- read.table(file_names[i], header = FALSE)
        colnames(data) <-c('chr','pos', 'N', 'X')

        # Store the data frame in the list
        data_list[[i]] <- data
    }

    BSobj <- makeBSseqData(data_list, sample_names)

    return(BSobj)
}

by_group_comparison <- function(bs_object, group_1_names, group_2_names, n_threads, output_prefix) {
    dmlTest <- DMLtest(bs_object, group1=group_1_names, group2=group_2_names, ncores=n_threads, smoothing=TRUE, smoothing.span=500)
    dmls <- callDML(dmlTest, delta=0.1, p.threshold=1e-5)
    dmrs <- callDMR(dmlTest, delta=0.1, p.threshold=1e-5)
    write.table(dmrs, file=paste(output_prefix, "DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(dmls, file=paste(output_prefix, "DML.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
}

create_hypothesis_table <- function(case_names, member) {
  df <- data.frame(case = case_names, sample = member)
  return(df)
}

use_hypothesis_testing <- function(bs_object, design, output_prefix) {

    design[] <- lapply(design, factor)

    DMLfit <- DMLfit.multiFactor(BSobj, design=design, formula=~case+origin+case:origin)

    DMLtest.hap <- DMLtest.multiFactor(DMLfit, coef=2)
    DMLtest.case <- DMLtest.multiFactor(DMLfit, coef=3)
    DMLtest.both <- DMLtest.multiFactor(DMLfit, coef=4)

    DMLtest.case <- DMLtest.multiFactor(DMLfit, coef=1)

    hap.ix <- sort(DMLtest.hap[,"pvals"], index.return=TRUE)$ix
    case.ix <- sort(DMLtest.case[,"pvals"], index.return=TRUE)$ix
    both.ix <- sort(DMLtest.both[,"pvals"], index.return=TRUE)$ix

    call.hap.dmr <- callDMR(DMLtest.hap, p.threshold=0.05)
    call.case.dmr <- callDMR(DMLtest.case, p.threshold=1e-5)
    call.both.dmr <- callDMR(DMLtest.both, p.threshold=0.05)

    write.table(call.hap.dmr, file=paste(output_prefix, "hap_DMR.tsv", sep="_"), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(call.case.dmr, file=paste(output_prefix, "case_DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(call.both.dmr, file=paste(output_prefix, "hap+case_DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
}

####### -------------- Analysis -------------- #######
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

file.names <- snakemake@input[["file_names"]]

print(file.names)

desired.analysis <- snakemake@wildcards[["analysis_type"]]
output.prefix <- snakemake@params[["output_prefix"]]
n.threads <- snakemake@threads

param_dict <- snakemake@params[[1]]

print(param_dict)

BSobj <- process_files(file.names, param_dict$sample_names)

if (desired.analysis == "two_group") {
    by_group_comparison(bs_object = BSobj, group_1_names = param_dict$groupA, group_2_names = param_dict$groupB, n_threads = n.threads, output_prefix=output.prefix)
} else if (desired.analysis == "model_based") {
    design.table <- create_hypothesis_table(case_names = param_dict$case, parental_origin = param_dict$parental_origin)
    design.table <- create_hypothesis_table(case_names = param_dict$case, member = param_dict$member)
    print(design.table)
    use_hypothesis_testing(bs_object = BSobj, design = design.table, output_prefix = output.prefix)
} else {
    warning(paste(desired.analysis, "has unsupported argument.", sep=" "))
}