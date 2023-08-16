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

test_model_one <- function(bs_object, design, output_prefix, which_base) {

    design[] <- lapply(design, factor)

    # Define the factor order
    if (which_base == "sibling") {
        design$case <- factor(design$case, levels=c("sibling", "mother", "asd"))
        coef1_name="casemother"
    else if (which_base == "mother") {
        design$case <- factor(design$case, levels=c("mother", "sibling", "asd"))
        coef1_name="casesibling"
    } else {
        warning(paste(which_base, "is an unsupported argument.", sep=" "))
        stop("Halting the code due to unsupported argument.")
    }

    DMLfit <- DMLfit.multiFactor(BSobj, design=design, formula=~case+family.id)

    print(DMLfit$X)
    cat("\n")

    DMLtest.case <- DMLtest.multiFactor(DMLfit, term="case")
    DMLtest.familyid <- DMLtest.multiFactor(DMLfit, term="family.id")
    DMLtest.case.coef1 <- DMLtest.multiFactor(DMLfit, coef=2)
    DMLtest.case.coef2 <- DMLtest.multiFactor(DMLfit, coef=3)

    case.ix <- sort(DMLtest.case[,"pvals"], index.return=TRUE)$ix
    familyid.ix <- sort(DMLtest.familyid[,"pvals"], index.return=TRUE)$ix
    case.coef1.ix <- sort(DMLtest.case.coef1[,"pvals"], index.return=TRUE)$ix
    case.coef2.ix <- sort(DMLtest.case.coef2[,"pvals"], index.return=TRUE)$ix

    write.table(DMLtest.case[case.ix,], file=paste(output_prefix, "case_DML.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(DMLtest.familyid[familyid.ix,], file=paste(output_prefix, "familyid_DML.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(DMLtest.case.coef1[case.coef1.ix,], file=paste(output_prefix, coef1_name, "DML.tsv", sep="_"), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(DMLtest.case.coef2[case.coef2.ix,], file=paste(output_prefix, "caseasd_DML.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)

    call.case.dmr <- callDMR(DMLtest.case, p.threshold=0.001)
    call.familyid.dmr <- callDMR(DMLtest.familyid, p.threshold=0.001)
    call.case.coef1.dmr <- callDMR(DMLtest.case.coef1, p.threshold=0.001)
    call.case.coef2.dmr <- callDMR(DMLtest.case.coef2, p.threshold=0.001)

    write.table(call.case.dmr, file=paste(output_prefix, "case_DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(call.familyid.dmr, file=paste(output_prefix, "familyid_DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(call.case.coef1.dmr, file=paste(output_prefix, coef1_name, "DMR.tsv", sep="_"), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(call.case.coef2.dmr, file=paste(output_prefix, "caseasd_DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
}

test_model_two <- function(bs_object, design, output_prefix) {

    design[] <- lapply(design, factor)

    DMLfit <- DMLfit.multiFactor(BSobj, design=design, formula=~case+allele+family.id)

    print(DMLfit$X)
    cat("\n")

    DMLtest.case <- DMLtest.multiFactor(DMLfit, term="case")
    DMLtest.familyid <- DMLtest.multiFactor(DMLfit, term="family.id")
    DMLtest.allele <- DMLtest.multiFactor(DMLfit, term="allele")

    case.ix <- sort(DMLtest.case[,"pvals"], index.return=TRUE)$ix
    familyid.ix <- sort(DMLtest.familyid[,"pvals"], index.return=TRUE)$ix
    allele.ix <- sort(DMLtest.allele[,"pvals"], index.return=TRUE)$ix

    write.table(DMLtest.case[case.ix,], file=paste(output_prefix, "case_DML.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(DMLtest.familyid[familyid.ix,], file=paste(output_prefix, "familyid_DML.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(DMLtest.allele[allele.ix,], file=paste(output_prefix, "allele_DML.tsv", sep="_"), sep="\t", quote=FALSE, row.names=FALSE)

    call.case.dmr <- callDMR(DMLtest.case, p.threshold=0.001)
    call.familyid.dmr <- callDMR(DMLtest.familyid, p.threshold=0.001)
    call.allele.dmr <- callDMR(DMLtest.allele, p.threshold=0.001)

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
    cat("\n")
    test_model_one(bs_object = BSobj, design = design.table, output_prefix = output.prefix, which_base=snakemake@wildcards[["which_base"]])
} else if (param_dict$model == "model_two") {
    design.table <- create_model_two_table(case_names = param_dict$case, allele_names= param_dict$allele_names, family_names = param_dict$family_names)
    print(design.table)
    cat("\n")
    test_model_two(bs_object = BSobj, design = design.table, output_prefix = output.prefix)
} else {
    warning(paste(param_dict$model, "has unsupported argument.", sep=" "))
    stop("Halting the code due to unsupported argument.")
}
