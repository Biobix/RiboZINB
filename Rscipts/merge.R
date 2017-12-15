############
#	R script merge zero inflated results
#	Outputs a file with best ISOFORM
############


#------ FUNCTIONS ------#
merge_files <- function(path,pattern){

	file_list <- list.files(path=path, pattern=paste(pattern,"$",sep=""))
 
	for (file in file_list) { 

	  	if (!exists("dataset")){ # if the merged dataset doesn't exist, create it
			dataset <- read.table(file=paste(path,file, sep="/"), header=TRUE, sep="\t", quote = "")
	  	} else {	# if the merged dataset does exist, append to it
			temp_dataset <-read.table(file=paste(path,file, sep="/"), header=TRUE, sep="\t", quote = "")
            #if(!exists("temp_dataset")){print(file)}
			dataset<-rbind(dataset, temp_dataset)
			rm(temp_dataset)
	 	}
	}

	return(dataset)
}


merge_files_2 <- function(path,pattern){

	file_list <- list.files(path=path, pattern=paste(pattern,"$",sep=""))

	for (file in file_list) {
	  	if (!exists("dataset")){ # if the merged dataset doesn't exist, create it
			dataset <- read.table(file=paste(path,file, sep="/"), header=F, sep="\t", quote = "")
	  	} else { # if the merged dataset does exist, append to it
			temp_dataset <-read.table(file=paste(path,file, sep="/"), header=F, sep="\t", quote = "")
			dataset<-rbind(dataset, temp_dataset)
			rm(temp_dataset)
	 	}
	}

	return(dataset)
}

#----------------------#


args <- commandArgs(TRUE)
dir <- as.character(args[1])
result_file <- as.character(args[2])
type <- as.character(args[3])


# positive set

if (type == "result") {

	transcript <- merge_files(dir,"_result")
	write.table(transcript, file=result_file, sep = "\t",col.names=T,row.names=F, quote = FALSE)

} else if (type == "perm")  {
	transcript.neg <- merge_files_2(dir,"_perm")
	write.table(transcript.neg, file=result_file, sep = "\t",col.names=F,row.names=F, quote = FALSE)
}




