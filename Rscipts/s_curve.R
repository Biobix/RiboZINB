############
#	R script to perform non parametric test
#	takes as input the directory of transcript files
#	Outputs a file with best ISOFORM
############


#----- Input Arguments -----#
args <- commandArgs(TRUE)
transcript_file <- as.character(args[1])
thresholds_file <- as.character(args[2])
scurve_file <- as.character(args[3])

# transcript data file
transcripts <- read.table(file=transcript_file,sep='\t',h=T)
transcripts <- transcripts[which(transcripts$rpkm_cds > 0),]
transcripts$log_rpkm <- log(transcripts$rpkm_cds)


# Implement S-Curve model on longest CCDS transcripts
transcripts.can <- transcripts[which(transcripts$CCDS == 'Y'),]
transcripts.can <- transcripts.can[order(transcripts.can$length_cds, transcripts.can$length_tr, decreasing=TRUE),] 


# USe the cononical/longest CCDS to estimate the threshold
CCDS <- character()

# if no CCDS information is available use the longest CCDS
if (nrow(transcripts.can) < 50) {
    all_genes <- as.vector(unique(transcripts$gene))
    for (i in 1:length(all_genes)) {
	    gene_tr <- transcripts[which(transcripts$gene == all_genes[i]),]
        CCDS[i] <- as.vector(gene_tr[order(gene_tr$length_cds, gene_tr$length_tr, decreasing=TRUE),]$transcript)[1]
    }
} else {
    all_genes <- as.vector(unique(transcripts.can$gene))
    for (i in 1:length(all_genes)) {
	    CCDS[i] <- as.vector(transcripts.can[which(transcripts.can$gene == all_genes[i]),]$transcript)[1]
    }
}

transcripts.scurve <- transcripts[which(transcripts$transcript %in% CCDS),]

# fit parameter logistic Sigmoid curve model
# A lower horizontal asymptote (minimum)
# B upper horizontal asymptote (maximum)
# xmid inflation point
# scal slope 

S.Model <- nls(coverage_cds ~ SSfpl(log_rpkm, A, B, xmid, scal), data=transcripts.scurve)

params <- coef(S.Model)
A <- params[1]
B <- params[2]
xmid <- params[3]
scal <- params[4]


# Bend point formula 2 [4 parameter logistic regression]
k <- 4.6804985798829056
Y_bend_lower = ((A - B)/(1 + (1/k))) + B
X_bend_lower = xmid*(((A - Y_bend_lower)/(Y_bend_lower - B))^(1/scal))

Y_bend_higher = ((A - B)/(1 + k)) + B
X_bend_higher = xmid*(((A - Y_bend_higher)/(Y_bend_higher - B))^(1/scal))

MINRPKM <- round(X_bend_lower, digits=4)
MINCOV <- predict(S.Model, newdata = data.frame(log_rpkm = X_bend_lower))[1]
MINCOV <- round(MINCOV, digits=4)

ax <- 2.5
pdf(file=scurve_file, width=10, useDingbats=FALSE)
mar.default = c(5, 4, 4, 2) + 0.1
par(mar=mar.default + c(0,3,0,0)) 
plot(transcripts.scurve$log_rpkm,transcripts.scurve$coverage_cds, main=paste(paste("Coverage",MINCOV, sep=" = "),paste("log RPKM",MINRPKM, sep=" = "), sep=",  "), ylab="RPF Coverage", xlab="Read Density [log RPKM]", cex.lab = ax, cex.axis=ax )
points(transcripts.scurve$log_rpkm,predict(S.Model), col='red', cex=1)
abline(h=MINCOV, col="blue", lwd=1.5)
abline(v=MINRPKM, col="blue", lwd=1.5)
invisible(dev.off())

MINRPKM <- exp(MINRPKM)
thresholds <- data.frame(coverage_cds = MINCOV, rpkm_cds = MINRPKM)

cat("CDS coverage threshold ", MINCOV, sep="", "\n")
cat("CDS read density threshold ", MINRPKM, sep="", "\n")

INFCOV <- predict(S.Model, newdata = data.frame(log_rpkm = xmid))[1]
cat("S curve inflation point (RPKM) ", xmid, sep="", "\n")
cat("S curve inflation point (coverage) ", INFCOV, sep="", "\n")


write.table(thresholds, file=thresholds_file, sep = "\t",col.names=T,row.names=F, quote = FALSE)


