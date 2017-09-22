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
all_genes <- as.vector(unique(transcripts.can$gene))
CCDS <- character()
for (i in 1:length(all_genes)) {
	CCDS[i] <- as.vector(transcripts.can[which(transcripts.can$gene == all_genes[i]),]$transcript)[1]
}

transcripts.scurve <- transcripts.can[which(transcripts.can$transcript %in% CCDS),]

# check if there are enough data points to estimant  s curve thresholds.
if (nrow(transcripts.scurve) <= 50) {

    cat("Not enough data points to generate S curve estimated thresholds default values will be used ", "", sep="", "\n")

    MINCOV_default <- 0.1
    MINRPKM_default <- 3
    cat("CDS coverage threshold ", MINCOV_default, sep="", "\n")
    cat("CDS read density threshold ", MINRPKM_default, sep="", "\n")

    thresholds_default <- data.frame(coverage_cds = MINCOV_default, rpkm_cds = MINRPKM_default)
    write.table(thresholds_default, file=thresholds_file, sep = "\t",col.names=T,row.names=F, quote = FALSE)
    stop()
}


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

pdf(file=scurve_file, useDingbats=FALSE)
plot(transcripts.scurve$log_rpkm,transcripts.scurve$coverage_cds, main=paste(paste("Coverage",MINCOV, sep=" = "),paste("log RPKM",MINRPKM, sep=" = "), sep=",  "), ylab="RPF Coverage", xlab="Read Density [log RPKM]" )
points(transcripts.scurve$log_rpkm,predict(S.Model), col='red')
abline(h=MINCOV, col="blue")
abline(v=MINRPKM, col="blue")
invisible(dev.off())

MINRPKM <- exp(MINRPKM)
thresholds <- data.frame(coverage_cds = MINCOV, rpkm_cds = MINRPKM)

cat("CDS coverage threshold ", MINCOV, sep="", "\n")
cat("CDS read density threshold ", MINRPKM, sep="", "\n")

INFCOV <- predict(S.Model, newdata = data.frame(log_rpkm = xmid))[1]
cat("S curve inflation point (RPKM) ", xmid, sep="", "\n")
cat("S curve inflation point (coverage) ", INFCOV, sep="", "\n")


write.table(thresholds, file=thresholds_file, sep = "\t",col.names=T,row.names=F, quote = FALSE)


