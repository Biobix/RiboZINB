############
#	R script to perform non parametric test
#	takes as input the directory of transcript files
#	Outputs a file with best ISOFORM
############

suppressMessages(library(ggplot2))
suppressMessages(library(DMwR))

args <- commandArgs(TRUE)
fdr <- as.numeric(args[1])
cutoff <- as.numeric(args[2])
default <- as.character(args[3])
transcript_file <- as.character(args[4])
prefix <- as.character(args[5])


###---- Data ----###
transcripts <- read.table(file=transcript_file,sep='\t',h=T, quote="")


# check if cutoff is estimated
if (default == 'y'){

    transcripts.FDR <- transcripts[which(transcripts$RPSS <= cutoff),]

} else if (default == 'n') {

    RPSS_pos <- na.omit(unique(sort(transcripts$RPSS, decreasing = FALSE)))
    RPSS_neg <- na.omit(unique(sort(transcripts$RPSS_neg, decreasing = FALSE)))
    RPSS <- RPSS_pos
    fdr_RPSS <- vector()
    for (i in 1:length(RPSS)) {
        RPSS.sel <- RPSS[i]
	    A <- length(RPSS_pos[RPSS_pos <= RPSS.sel])
	    B <- length(RPSS_neg[RPSS_neg <= RPSS.sel])
	    fdr_RPSS[i] <- B/(A + B)
    }

    FDR_RPSS <- data.frame(RPSS=RPSS,fdr=fdr_RPSS)
    cutoff <- max(FDR_RPSS[which(FDR_RPSS$fdr <= fdr),]$RPSS)

    cutoff <- round(cutoff + 0.005, digits=2)    # round-up to the nearest 3rd decimal

    transcripts.FDR <- transcripts[which(transcripts$RPSS <= cutoff),]

    # plot positive and negative RPSS distribution
    xlim <- 5
    pdf(file=paste(prefix,"distributions.pdf",sep="_"), width=10, useDingbats=FALSE)
    ax <- 30
    dta.RPSS <- data.frame(RPSS = c(RPSS_pos,RPSS_neg), Distribution=c(rep("Observed RPSS",length(RPSS_pos)),rep("Empirical RPSS",length(RPSS_neg)))) 
    print(ggplot(dta.RPSS, aes(x = RPSS, fill = Distribution)) + geom_density(alpha = 0.5, na.rm=T)  + theme_bw(base_size=ax) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + xlim(0, xlim)  + theme(legend.position = c(0.5, 0.8), text = element_text(size=ax), legend.key.size = unit(1, "cm"),legend.title.align=1, axis.text.y = element_text(angle = 90, vjust=0, hjust=0.5), axis.title = element_text(colour = "black")))
    invisible(dev.off())
}

transcripts.FDR <- transcripts[which(transcripts$RPSS <= cutoff),]

cat("\nTotal number of genes expressed", length(unique(transcripts$gene)), sep=" ", "\n")
cat("Total number of transcripts", length(unique(transcripts$transcript)), sep=" ", "\n")

cat("\nThreshold RPSS", cutoff, sep=" ", "\n")
cat("Number of genes below threshold", length(unique(transcripts.FDR$gene)), sep=" ", "\n")
cat("Number of trancript below", length(unique(transcripts.FDR$transcript)), sep=" ", "\n")

selected_genes <- as.vector(unique(transcripts.FDR$gene))
genes_single <- character()
genes_multi <- character()
p <- 1
t <- 1
for (i in 1:length(selected_genes)) {
	gene_tr <- nrow(transcripts.FDR[which(transcripts.FDR$gene == selected_genes[i]),])
	if (gene_tr == 1) {
		genes_single[p] <- selected_genes[i]
		p <- p + 1
	} else if (gene_tr > 1) {
		genes_multi[t] <- selected_genes[i]
		t <- t + 1
	}
}

cat("Number of genes with single isoform ", length(genes_single), sep=" ", "\n")
cat("Number of genes with multiple isoform ", length(genes_multi), sep="", "\n")

# add best scoring isoform for all those below cutoff
gene_no_tr_below_cutoff <- setdiff(transcripts$gene,transcripts.FDR$gene)
cat("Number of genes without any isoforms below threshold ", length(gene_no_tr_below_cutoff), sep="", "\n")


# write output to file
expressed_file <- paste(prefix,"expressed_isoforms.txt",sep="_")
write.table(transcripts.FDR, file=expressed_file, sep = "\t",col.names=T,row.names=F, quote = FALSE)


