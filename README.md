# RiboZINB
RiboZinB: A tool to identifying actively translated isoform(s) from ribosome profiling data. The case of the Zero-inflated model.

The software was created at VIB-UGent Center for Medical Biotechnology and Lab of Bioinformatics and Computational Genomics (Biobix), University of Gent, belgium.


# Requirements

The RiboZINB software is primarily written in perl and R.

	Perl packages:
	--------------
		- Perl + v.5.10+
		- Getopt::Long
		- Cwd	
		- BioPerl
		- POSIX
		- Parallel::ForkManager
		- File::Spec
		
	R packages:
	--------------
		- ggplot2
		- pscl
		- VGAM

	
# Install

RiboZINB does not rerquire any special installation requirements. To install RiboZINB, clone/download the github repository into a directory of your choice.  

USAGE: perl ./RiboZINB.pl -p positional_data_file -g gtf_file -s script_dir -r total_mappable_reads  

Mandatory input parameters  
	-p	RIBO-seq positional data in comma delimited format [chromosome,strand,start,count]  
	-g	annoation in gtf format  
	-s	script directory  
	-r	total number of mappable reads  

Optional input parameters  
	-e	experiment name  
	-w	work directory  
	-m	minimum reads count  
	-d	delete work directory if it exist  
	-t	number of processors to use  
	-n	number of iterations when generating negative set  
	-f	FDR threshold  
    -v	Default RPSS threshold [v=0.1]  
	-dt flag to determine if we use default RPSS threshold [0.1] or estimate threshold by permutation test.  
    -a	multiplier to adjust for noise  


Output files  
	*_gene_list.txt	this file contains the list of all protein coding genes evaluated in the analysis  
	*_all_transcript.txt File consisting of all protein coding transcripts from the genes in *_gene_list.txt  
	*_transcript_zeroinfl.txt	file consisting of transcripts from _all_transcript.txt that pass the S-curve thresholds alongside their RiboZinB scores and P-values  
	*_Ribo_expressed_isoform.txt consist of the list of transcripts/isoforms predicted by the RiboZinB algorithm i.e. transcripts that fall below the FDR threshold  
	*_Ribo_distribution.pdf	plot of the empirical and observed distributions of the RPSS scores  
	*_scurve.pdf	plot of the fitted s-curve model overlay of on the longest CCDS (canonical isoform) and the log RPKM and CDS coverage thresholds  

	
# More information  
For more information about RiboZINB contact: elvis.ndah@gmail.com or Gerben.Menschaert@UGent.be

