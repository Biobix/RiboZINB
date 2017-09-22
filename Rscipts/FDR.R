############
#	R script to perform non parametric test
#	takes as input the directory of transcript files
#	Outputs a file with best ISOFORM
############

suppressMessages(library(ggplot2))

args <- commandArgs(TRUE)
fdr <- as.numeric(args[1])
cutoff <- as.numeric(args[2])
default <- as.character(args[3])
transcript_file <- as.character(args[4])
prefix <- as.character(args[5])
fdr_type <- as.character(args[6]) 


###---- Data ----###
transcripts <- read.table(file=transcript_file,sep='\t',h=T, quote="")

# check if cutoff is estimated
if (default == 'y'){

    transcripts.FDR <- transcripts[which(transcripts$score <= cutoff),]

} else if (default == 'n') {

    scores_pos <- na.omit(sort(transcripts$score, decreasing = FALSE))
    scores_neg <- na.omit(sort(transcripts$score_neg, decreasing = FALSE))

    scores <- na.omit(sort(transcripts$score, decreasing = FALSE))
    fdr_score <- vector()
    for (i in 1:length(scores)) {
        score.sel <- scores[i]
	    A <- length(scores_pos[scores_pos <= score.sel])
	    B <- length(scores_neg[scores_neg <= score.sel])
	    fdr_score[i] <- B/(A + B)
    }

    FDR_score <- data.frame(score=scores,fdr=fdr_score)
    cutoff <- max(FDR_score[which(FDR_score$fdr <= fdr),]$score)

    transcripts.FDR <- transcripts[which(transcripts$score <= cutoff),]

    # plot positive and negative score distribution
    xlim <- 5
    ppi <- 300
    png(paste(prefix,"distributions.png",sep="_"), width=6*ppi, height=6*ppi, res=ppi)
    #pdf(file=paste(prefix,"distributions.pdf",sep="_"), width=12, useDingbats=FALSE)

    dta.score <- data.frame(scores = c(scores_pos,scores_neg), Distribution=c(rep("Observed",length(scores_pos)),rep("Negative",length(scores_neg)))) 
    print(ggplot(dta.score, aes(x = scores, fill = Distribution)) + geom_density(alpha = 0.5, na.rm=T)  + theme_bw(base_size=18) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + xlim(0, xlim)  + theme(legend.position = c(0.5, 0.8), text = element_text(size=20), legend.key.size = unit(0.5, "cm"),legend.title.align=0.2))


}


transcripts.FDR <- transcripts[which(transcripts$score <= cutoff),]
transcripts.FDR$isoform_above_cutoff <- "Y"


cat("Total number of genes expressed", length(unique(transcripts$gene)), sep=" ", "\n")
cat("Total number of transcripts", length(unique(transcripts$transcript)), sep=" ", "\n")

cat("\nThreshold score", cutoff, sep=" ", "\n")
cat("Number of genes above threshold", length(unique(transcripts.FDR$gene)), sep=" ", "\n")
cat("Number of  trancript called", length(unique(transcripts.FDR$transcript)), sep=" ", "\n")

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
	} else {
		genes_multi[t] <- selected_genes[i]
		t <- t + 1
	}
}

cat("Number of genes with single isoform ", length(genes_single), sep=" ", "\n")
cat("Number of genes with multiple isoform ", length(genes_multi), sep="", "\n")
cat("","",sep="","\n")

# add best scoring isoform for all those below cutoff
cat("Select best scoring for all genes that did not cross the thresholds ", cutoff, sep="", "\n")
gene_no_tr_above_cutoff <- setdiff(transcripts$gene,transcripts.FDR$gene)


if (length(gene_no_tr_above_cutoff) > 0) {

    transcripts.low <- transcripts[which(transcripts$gene %in% gene_no_tr_above_cutoff),]
    transcripts.low$isoform_above_cutoff <- "N"

    best_isoform <- character()
    selected_genes <- as.vector(unique(transcripts.low$gene))

    t <- 1
    for (i in 1:length(selected_genes)) {

	    tr <- transcripts.low[which(transcripts.low$gene == selected_genes[i]),]

	    tr_id <- as.vector(tr[which(tr$score == min(tr$score)),]$transcript)
        if (length(tr_id) > 1) {
            for (j in 1:length(tr_id)) {
                best_isoform[t] <- tr_id[j]
                t <- t + i
            }
        } else {
	        best_isoform[t] <- tr_id[1]
            t <- t + 1
        }
    }

    transcripts.low.sel <- transcripts.low[which(transcripts.low$transcript %in% best_isoform),]
    transcripts.FDR <- rbind(transcripts.low.sel,transcripts.FDR)
}


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
	} else {
		genes_multi[t] <- selected_genes[i]
		t <- t + 1
	}
}

cat("Number of genes called", length(unique(transcripts.FDR$gene)), sep=" ", "\n")
cat("Number of  trancript called", length(unique(transcripts.FDR$transcript)), sep=" ", "\n")

cat("","",sep="","\n")
cat("Number of genes with single isoform ", length(genes_single), sep=" ", "\n")
cat("Number of genes with multiple isoform ", length(genes_multi), sep="", "\n")


composite_merge_file <- paste(prefix,"merged_expressed_isoforms.txt",sep="_")
write.table(transcripts.FDR, file=composite_merge_file, sep = "\t",col.names=T,row.names=F, quote = FALSE)


#### Best scoring isoforms
selected_genes <- as.vector(unique(transcripts$gene))
best_scoring_tr <- vector()
t <- 1
for (i in 1:length(selected_genes)) {
	gene_tr <- transcripts[which(transcripts$gene == selected_genes[i]),]

    tr_list <- as.vector(gene_tr[which(gene_tr$score == min(gene_tr$score)),]$transcript)

    for (j in 1:length(tr_list)) {
        best_scoring_tr[t] <- tr_list[j]
        t <- t + 1
    }
}


print(tr_list)

transcripts.best_isoform <- transcripts[which(transcripts$transcript %in% best_scoring_tr),]
best_scoring_isoform <- paste(prefix,"best_scoring_isoform.txt",sep="_")
write.table(transcripts.best_isoform, file=best_scoring_isoform, sep = "\t",col.names=T,row.names=F, quote = FALSE)

cat("Number of unique CCDS_IDs in best isoform list ", length(unique(transcripts.best_isoform$CCDS_ID)), sep="", "\n")




