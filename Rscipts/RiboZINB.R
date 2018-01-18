############
#	R script to perform non parametric test
#	takes as input the directory of transcript files
#	Outputs a file with best ISOFORM
############
options(warn=-1)
suppressMessages(library(pscl))
suppressMessages(library(VGAM))


#------------ FUNCTIONS ------------#
zero_inflated <- function(frame){
	ctrl <- zeroinfl.control(method = "CG", reltol=-Inf)
	zinb <- zeroinfl(count ~., data = frame, link = "logit", dist = "negbin",control = ctrl)
	return(zinb)
}

ZINB <- function(frame){
  return(tryCatch(zero_inflated(frame),
     error = function(e){}))
}

# Fisher yates shuffling algorithm
fisher_yates <- function(counts) {
    n <- length(counts)
    i <- n
    while (i > 0) {
        j <- floor(runif(1)*i+1)
        if (i != j) {
            tmp <- counts[i]
            counts[i] <- counts[j]
            counts[j] <- tmp
        }
        i <- i - 1
    }
    return(counts)
}


#----- Input Arguments -----#
args <- commandArgs(TRUE)
dir <- as.character(args[1])
transcript_file <- as.character(args[2])
result_file <- as.character(args[3])
thresholds_file <- as.character(args[4])
no_of_samples <- as.numeric(args[5])
alpha <- as.numeric(args[6])


pause_threshold <- 0.6  # pause threshold proportion

# path to data files
files <- list.files(path=dir)
path <- paste(dir, "/", sep="")

# transcript data file
transcript.dta <- read.table(file=transcript_file,sep='\t',h=T, quote="")
thresholds <- read.table(file=thresholds_file,sep='\t',h=T)

# keep only transcript with some evidence of ribosome occupancy
transcript.dta <- transcript.dta[which(transcript.dta$rpkm_cds > 0),]
transcript.dta <- transcript.dta[which(transcript.dta$coverage_cds > 0),]

# get genes with at least one transcript above threshold and perform predcitions on all transcripts in the gene model
transcript.threshold <- transcript.dta[which(transcript.dta$coverage_cds >= thresholds$coverage_cds),]
transcript.threshold <- transcript.threshold[which(transcript.threshold$rpkm_cds >= thresholds$rpkm_cds),]

transcript.dta <- transcript.dta[which(transcript.dta$gene %in% transcript.threshold$gene),]

# new variables to update
transcript.dta$rate_ratio <- vector(length = nrow(transcript.dta))
transcript.dta$odds_ratio <- vector(length = nrow(transcript.dta))
transcript.dta$RPSS <- vector(length = nrow(transcript.dta))
transcript.dta$RPSS_neg <- vector(length = nrow(transcript.dta))
transcript.dta$LCI_count <- vector(length = nrow(transcript.dta))
transcript.dta$UCI_count <- vector(length = nrow(transcript.dta))
transcript.dta$LCI_zero <- vector(length = nrow(transcript.dta))
transcript.dta$UCI_zero <- vector(length = nrow(transcript.dta))
transcript.dta$pvalues_count <- vector(length = nrow(transcript.dta))
transcript.dta$pvalues_zero <- vector(length = nrow(transcript.dta))


permutations <- list()

for (t in 1:nrow(transcript.dta)) {

	gene.file <- paste(transcript.dta[t,]$gene,".txt", sep="")
	gene <- read.table(file=paste(path,gene.file, sep=""),sep='\t',h=T, quote="")
	gene.f <- data.frame(count=gene[,1], group='1')

	j <- which(transcript.dta[t,]$transcript == colnames(gene)) # get column number of transcript at row t
	tr.df <- data.frame(count=gene[,j], group='0')
    count_values <- gene[,1]

    # adjust for spikes/pause sites by replacing positions with reads > 70% with average
    #average_reads <- round(0.5 + sum(count_values)/length(count_values[count_values > 0]), digits=0)
    #pause_site_threshold <- pause_threshold*sum(count_values)       # threshold for pause sites

    # replace pause sites with average reads
    #gene.f$count[gene.f$count >= pause_site_threshold] <- average_reads
    #tr.df$count[tr.df$count >= pause_site_threshold] <- average_reads

    # inflate zero if needed uninflated zeros occurs mostly within histon genes
    proportion_of_count <- length(count_values[count_values > 0])/length(count_values)
    if (proportion_of_count > 0.5) {
        adjusted_zeros <- round(proportion_of_count*nrow(gene), digits=0)
	    gene.f <- data.frame(count=c(gene[,1],rep(0,adjusted_zeros)), group='1')
	    tr.df <- data.frame(count=c(gene[,j],rep(0,adjusted_zeros)), group='0')
    }

    # create combined data frame
	frame <- rbind(gene.f,tr.df)
	frame$group <- as.factor(frame$group)
	frame$count <- round(frame$count, 0)

	zinb <- ZINB(frame)

	if (is.null(zinb) == 'FALSE') {

        #estimates
        rate_ratio  <- exp(coef(zinb)[2])
        odds_ratio <- exp(coef(zinb)[4])

        transcript.dta[t,]$rate_ratio <- rate_ratio
        transcript.dta[t,]$odds_ratio <- odds_ratio

		# RPSS
        x <- c(1 - rate_ratio, 1 - odds_ratio)
		transcript.dta[t,]$RPSS <- sqrt(sum(x^2))

		# P value
		transcript.dta[t,]$pvalues_count <- summary(zinb)$coefficients$count[2,4]
		transcript.dta[t,]$pvalues_zero <- summary(zinb)$coefficients$zero[2,4]

		# Confidence interval
		transcript.dta[t,]$LCI_count <- exp(confint(zinb)[2,1])
		transcript.dta[t,]$UCI_count <- exp(confint(zinb)[2,2])
		transcript.dta[t,]$LCI_zero <- exp(confint(zinb)[4,1])
		transcript.dta[t,]$UCI_zero <- exp(confint(zinb)[4,2])

        # Permutation test
	    RPSS_neg.samp <- vector()
	    for (i in 1:no_of_samples) {

            # Simulated distribution
            set.seed(i)
            lambda <- predict(zinb, type = "count")
            p <- predict(zinb, type = "zero")

            mu <- alpha*mean(lambda)
            gene.r <- rzinegbin(nrow(gene.f), pstr0 = p[1:nrow(gene.f)], mu = mu, size = zinb$theta)

            gene_neg.df <- data.frame(count=gene.r + gene.f$count, group="1")
		    frame.r <- rbind(gene_neg.df,tr.df)
		    frame.r$group <- as.factor(frame.r$group)

		    zinb_neg <- ZINB(frame.r)
		    if (is.null(zinb_neg) == 'FALSE') {

			    rate_ratio_neg <- exp(coef(zinb_neg)[2])
			    odds_ratio_neg <- exp(coef(zinb_neg)[4])

                x.neg <- c(1 -rate_ratio_neg, 1 - odds_ratio_neg)
                RPSS_neg.samp[i] <- sqrt(sum(x.neg^2))

                rm(rate_ratio_neg)
                rm(odds_ratio_neg)
                rm(x.neg)
		    } else {
                RPSS_neg.samp[i] <- 'NA'
            }
	    }

        # write permutation estimates to file
        #permutations[[colnames(gene)[j]]] <- RPSS_neg.samp
		transcript.dta[t,]$RPSS_neg <- mean(RPSS_neg.samp, na.rm =T)
        
        rm(rate_ratio)
        rm(odds_ratio)

	} else {

		transcript.dta[t,]$RPSS  <- NA
		transcript.dta[t,]$pvalues_count <- NA
		transcript.dta[t,]$pvalues_zero <- NA
		transcript.dta[t,]$LCI_count <- NA
		transcript.dta[t,]$UCI_count <- NA
		transcript.dta[t,]$LCI_zero <- NA
		transcript.dta[t,]$UCI_zero <- NA
        transcript.dta[t,]$rate_ratio  <- NA
        transcript.dta[t,]$odds_ratio <- NA
        transcript.dta[t,]$RPSS_neg <- NA
	}

}

write.table(transcript.dta, file=paste(result_file, "_result", sep=""), sep = "\t",col.names=T,row.names=F, quote = FALSE)


# Write permutation list to file
#dataset.perm <- c()
#for (i in 1:length(permutations)) {
#    row <- c(names(permutations)[i], permutations[[i]])
#    dataset.perm <- rbind(dataset.perm, row)
#}
#dataset.perm <-as.data.frame(dataset.perm)
#colnames(dataset.perm) <- c("transcript", paste("perm_", c(1:no_of_samples),sep=""))
#write.table(dataset.perm, file=paste(result_file, "_perm", sep=""), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



