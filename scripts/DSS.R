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

create_model_one_table <- function(case_names, family_names) {
  df <- data.frame(case = case_names, family.id = family_names, stringsAsFactors = TRUE)
  return(df)
}

create_model_two_table <- function(case_names, allele_names, family_names) {
  df <- data.frame(case = case_names, allele = allele_names, family.id = family_names, stringsAsFactors = TRUE)
  return(df)
}

test_model_one <- function(bs_object, design, output_prefix) {

    design[] <- lapply(design, factor)

    DMLfit <- DMLfit.multiFactor(BSobj, design=design, formula=~case+family.id)

    print(DMLfit$X)

    DMLtest.case <- DMLtest.multiFactor(DMLfit, term="case")
    DMLtest.familyid <- DMLtest.multiFactor(DMLfit, term="family.id")
    DMLtest.case.coef1 <- DMLtest.multiFactor(DMLfit, coef=2)
    DMLtest.case.coef2 <- DMLtest.multiFactor(DMLfit, coef=3)


    case.ix <- sort(DMLtest.case[,"pvals"], index.return=TRUE)$ix
    familyid.ix <- sort(DMLtest.familyid[,"pvals"], index.return=TRUE)$ix
    case.coef1.ix <- sort(DMLtest.case.coef1[,"pvals"], index.return=TRUE)$ix
    case.coef2.ix <- sort(DMLtest.case.coef2[,"pvals"], index.return=TRUE)$ix

    call.case.dmr <- callDMR(DMLtest.case, p.threshold=0.05)
    call.familyid.dmr <- callDMR(DMLtest.familyid, p.threshold=0.05)
    call.case.coef1.dmr <- callDMR(DMLtest.case.coef1, p.threshold=0.05)
    call.case.coef2.dmr <- callDMR(DMLtest.case.coef2, p.threshold=0.05)

    write.table(call.case.dmr, file=paste(output_prefix, "case_DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(call.familyid.dmr, file=paste(output_prefix, "familyid_DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(call.case.coef1.dmr, file=paste(output_prefix, "case-coef1_DMR.tsv", sep="_"), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(call.case.coef2.dmr, file=paste(output_prefix, "case-coef2_DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
}

test_model_two <- function(bs_object, design, output_prefix) {

    design[] <- lapply(design, factor)

    DMLfit <- DMLfit.multiFactor(BSobj, design=design, formula=~case+allele+family.id)

    print(DMLfit$X)

    DMLtest.case <- DMLtest.multiFactor(DMLfit, term="case")
    DMLtest.familyid <- DMLtest.multiFactor(DMLfit, term="family.id")
    DMLtest.allele <- DMLtest.multiFactor(DMLfit, term="allele")

    case.ix <- sort(DMLtest.case[,"pvals"], index.return=TRUE)$ix
    familyid.ix <- sort(DMLtest.familyid[,"pvals"], index.return=TRUE)$ix
    allele.ix <- sort(DMLtest.allele[,"pvals"], index.return=TRUE)$ix

    call.case.dmr <- callDMR(DMLtest.case, p.threshold=0.05)
    call.familyid.dmr <- callDMR(DMLtest.familyid, p.threshold=0.05)
    call.allele.dmr <- callDMR(DMLtest.allele, p.threshold=0.05)

    write.table(call.case.dmr, file=paste(output_prefix, "case_DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(call.familyid.dmr, file=paste(output_prefix, "familyid_DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(call.allele.dmr, file=paste(output_prefix, "allele_DMR.tsv", sep="_"), sep="\t", quote=FALSE, row.names=FALSE)
}

####### -------------- Analysis -------------- #######
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

file.names <- snakemake@input[["file_names"]]

output.prefix <- snakemake@params[["output_prefix"]]
n.threads <- snakemake@threads

param_dict <- snakemake@params[[1]]

BSobj <- process_files(file.names, param_dict$sample_names)

if (param_dict$model == "model_one") {
    design.table <- create_model_one_table(case_names = param_dict$case, family_names = param_dict$family_names)
    print(design.table)
    test_model_one(bs_object = BSobj, design = design.table, output_prefix = output.prefix)
} else if (param_dict$model == "model_two") {
    design.table <- create_model_two_table(case_names = param_dict$case, allele_names= param_dict$allele_names, family_names = param_dict$family_names)
    print(design.table)
    test_model_two(bs_object = BSobj, design = design.table, output_prefix = output.prefix)
} else {
    warning(paste(param_dict$model, "has unsupported argument.", sep=" "))
}