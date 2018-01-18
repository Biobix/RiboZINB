#!/usr/bin/perl -w

use strict;
use warnings;
use diagnostics;

use DBI;
use Getopt::Long;
use POSIX;
use Parallel::ForkManager;
use File::Spec;
use POSIX qw(ceil);
use Cwd 'abs_path';


########
# USAGE: perl RiboZINB.pl -p positional_data_file -g gtf_file -e mESC -t tmp_folder -m minimum_reads_allowed -r total_mappable_reads
########


my $startRun = time();

my $pause_site_prop = 0.7;  #proportion of reads at one position to consider are pause site

my $default_score = "Y";# Use default score threshold [0.1] or estimate threshold by permutation test.
my $cutoff = 0.1;        # Default score threshold to determine expressed isoform(s)
my $MINCOUNT = 5;		# minimum reads count
my $force = "Y";		# Y to delete tmp folder if it exist
my $exp_name;			# experiment name
my $work_dir;			# working folder
my $ribo;				# RIBO-seq derived master SQLite file
my $gtf;				# GTF file
my $script_dir;		    # Directory where the script files are stored (defaults to the current directory)
my $mapped_total;		# total mappable reads
my $cores = 1;		    # number of threads
my $no_of_samples = 30;# number of iterations when generating negative set
my $fdr = 0.05;		    # FDR percentage
my $sigma = 1;


# Output files
my $transcript_info_file;
my $isoform;
my $scurve;
my $threshold;
my $transcript_threshold;


# Global variables hardcoded
my $space_len = 55;		# spaces to format output


##	GET command line arguments
GetOptions(

	'p=s'=> \$ribo,
	'g=s'=> \$gtf,
	'e=s'=> \$exp_name,
	'w=s'=> \$work_dir,
	'm=f'=> \$MINCOUNT,
	'r=i'=> \$mapped_total,
	'd=s'=> \$force,
	't=i'=> \$cores,
	'n=i'=> \$no_of_samples,
	's=s'=> \$script_dir,
	'f=f'=> \$fdr,
    'v=f'=> \$cutoff,
	'dt=s'=> \$default_score,
    'a=f'=>\$sigma,

    # Output files
    'tf=s'=>\$transcript_info_file,
    'pf=s'=>\$isoform,
    'sc=s'=>\$scurve,
    'th=s'=>\$threshold
);


#----------------------------------------------------
#			EXECUTION
#----------------------------------------------------

my %params = (
	'p'=> $ribo,
	'g'=> $gtf,
	's'=> $script_dir,
	'r'=> $mapped_total
);


#---- check if mandatory parameters are initialized
my @invalid = grep uninitialized_param($params{$_}), keys %params;	
die "Mandatory variable(s) not properly initialized: @invalid\n" if @invalid;


#---- Check working directory
my $file_dir;
my $tmp_dir;
($work_dir, $file_dir, $tmp_dir) = check_work_dir($work_dir);


#---- convert to uppercase
$default_score=uc($default_score);
$force=uc($force);


#----- Check if input varaibles are properly initialized
if (_isNumeric($cutoff) == 0 or (_isNumeric($cutoff) == 1 and ($cutoff <= 0 and $cutoff >= 1))) {
    print "INPUT ERROR: -v is not properly initialized!!!!!!!.\n";
    print "-v take as input a numeric value in the interval (0, 1).\n";
    exit(1);
}

if (_isNumeric($fdr) == 0 or (_isNumeric($fdr) == 1 and ($fdr <= 0 and $fdr >= 1))) {
    print "INPUT ERROR: -f is not properly initialized!!!!!!!.\n";
    print "-f take as input a numeric value in the interval (0, 1).\n";
    exit(1);
}

if (_isNumeric($MINCOUNT) == 0 or (_isNumeric($MINCOUNT) == 1 and $MINCOUNT < 1)) {
    print "INPUT ERROR: -m is not properly initialized!!!!!!!.\n";
    print "-m take as input a number greater than or equal to 1.\n";
    exit(1);
}

unless ($force ne "N" or $force ne "Y") {
    print "INPUT ERROR: -d is not properly initialized!!!!!!!.\n";
    print "-d take as input either N(n) or Y(y).\n";
    exit(1);
}

if (_isNumeric($mapped_total) == 0 or (_isNumeric($mapped_total) == 1 and $mapped_total <= 1)) {
    print "INPUT ERROR: -r is not properly initialized!!!!!!!.\n";
    print "-r take as input an integer greater than 1.\n";
    exit(1);
}

if (_isNumeric($cores) == 0 or (_isNumeric($cores) == 1 and $cores < 1)) {
    print "INPUT ERROR: -t is not properly initialized!!!!!!!.\n";
    print "-t take as input an integer greater than or equal to 1.\n";
    exit(1);
}

if (_isNumeric($no_of_samples) == 0 or (_isNumeric($no_of_samples) == 1 and $no_of_samples < 1)) {
    print "INPUT ERROR: -n is not properly initialized!!!!!!!.\n";
    print "-n take as input an integer greater than or equal to 1.\n";
    exit(1);
}

if (_isNumeric($no_of_samples) == 0 or (_isNumeric($no_of_samples) == 1 and $no_of_samples < 1)) {
    print "INPUT ERROR: -n is not properly initialized!!!!!!!.\n";
    print "-n take as input an integer greater than or equal to 1.\n";
    exit(1);
}

if (_isNumeric($sigma) == 0 or (_isNumeric($sigma) == 1 and $sigma <= 0)) {
    print "INPUT ERROR: -a is not properly initialized!!!!!!!.\n";
    print "-a take as input a numeric value greater than 0.\n";
    exit(1);
}

if ($script_dir) {
	if (-d $script_dir) {
		my $full_path = abs_path($script_dir);
		printf("%-".$space_len."s", "The following script directory is used ");
		print ":$full_path\n";

		unless (-e $script_dir."/RiboZINB.R") {
			print "Script 'RiboZINB.R' not found in script directory.\nEnsure the directory  '$script_dir' contains all required file (see readme).\n";
			exit(1);
		}

		unless (-e $script_dir."/merge.R") {
			print "Script 'merge.R' not found in script directory.\nEnsure the directory '$script_dir' contains all required file (see readme).\n";
			exit(1);
		}

		unless (-e $script_dir."/FDR.R") {
			print "Script 'FDR.R' not found in script directory.\nEnsure the script directory '$script_dir' contains all required file (see readme).\n";
			exit(1);
		}

		unless (-e $script_dir."/s_curve.R") {
			print "Script 's_curve.R' not found in script directory.\nEnsure the script directory '$script_dir' contains all required file (see readme).\n";
			exit(1);
		}
	} else { 
		print "Script directory '$script_dir' does NOT exist!\n";
		exit;
	}
} else {
	print "Please do not forget to pass the scripts directory!\n";
	exit;
}


#---- print input variables 
print "\n\n";
if ($exp_name) {
	printf("%-".$space_len."s", "Experiment name");
	print ":$exp_name\n";
    $exp_name = $exp_name."_";
} else {
	$exp_name = "riboZinB_pred_";
}

my $full_path = abs_path($work_dir);
printf("%-".$space_len."s", "Working directory",);
print ":$full_path\n";

if ($ribo) {
	if (-e $ribo) {
		my $full_path = abs_path($ribo);
		printf("%-".$space_len."s", "File containing RPF counts per genomic position");
		print ":$full_path\n";
	} else { 
		print "File '$ribo' does NOT exist!\n";
		exit;
	}
} else {
	print "Please do not forget to pass file with RPF RPF counts per genomic position!\n";
	exit;
}

if ($gtf) {
	if (-e $gtf) {
		my $full_path = abs_path($gtf);
		printf("%-".$space_len."s", "The GTF annotation file used");
		print ":$full_path\n";
	} else { 
		print "File '$gtf' does NOT exist!\n";
		exit;
	}
} else {
	print "Please do not forget to pass the GTF annotation file!\n";
	exit;
}

printf("%-".$space_len."s", "File containing RPF counts per genomic position");
print ":$full_path\n";


printf("%-".$space_len."s", "Total number of mappable reads",);
print ":$mapped_total\n";

printf("%-".$space_len."s", "Minimum read count allowed at the gene level",);
print ":$MINCOUNT\n";

printf("%-".$space_len."s", "Total number of cores used",);
print ":$cores\n";

if (uc($default_score) eq 'Y') {
    printf("%-".$space_len."s", "Threshold score used",);
    print ":$cutoff\n";
} else {
    printf("%-".$space_len."s", "False discovery rate used",);
    print ":$fdr\n";
}

printf("%-".$space_len."s", "Number of iterations to generate negative set ",);
print ":$no_of_samples\n";

printf("%-".$space_len."s", "Proportion of noise to generate negative set ",);
print ":$sigma\n";


#---- round non integers to integer values
$mapped_total = ceil($mapped_total);
$cores = ceil($cores);
$no_of_samples = ceil($no_of_samples);


#----- Prepare output files
my $gene_file = $work_dir."/".$exp_name."gene_list.txt";
$transcript_info_file = $work_dir."/".$exp_name."all_transcript.txt";
$isoform = $work_dir."/".$exp_name."transcript_zeroinfl.txt";
$scurve = $work_dir."/".$exp_name."scurve.pdf";
$threshold = $work_dir."/".$exp_name."thresholds.txt";
$transcript_threshold = $work_dir."/".$exp_name."transcript_threshold.txt";


#---- get all the genes and transcripts in the genome
print "\n\n";
my $msg = "\tExtracting annotation information from GTF file(s).";
print_localtime($msg);
my ($genes, $transcripts) = annotation_GTF($gtf);


#---- get read table into memory
$msg = "\tExtracting RPF positional data from file.";
print_localtime($msg);
my $reads_count_table = get_read_table($ribo);


#---- assemble transcripot information
$msg = "\tAssembling transcript information";
print_localtime($msg);
my ($valid_genes,$transcripts_valid) = Isoform_data($genes,$transcripts);


#---- write gene and transcript info to file
gene_info($valid_genes,$gene_file);
transcript_info($transcripts_valid,$transcript_info_file);


#---- use S curve to find cutoffs
$msg = "\tGenerating S curve thresholds.";
print_localtime($msg);
my $scurve_cmd = "Rscript $script_dir"."/s_curve.R $transcript_info_file $threshold $scurve";
print "\n$scurve_cmd\n";
system($scurve_cmd);


#---- split transcript into multiple files for multiprocessing
split_genes_file($transcripts_valid);


#---- Peroform zero inflated negative binomial regression
print "\n\n";
$msg = "\tPerforming Zeroinflated negative binomial analysis. This might take a while.";
print_localtime($msg);
ZINB_parallel();


#--- Estimate FDR and select expressed Isoforms
print "\n\n";
$msg = "\tIdentifing expressed isoform(s) at $fdr FDR threshold";
print_localtime($msg);


my $merge_cmd = "Rscript $script_dir"."/merge.R $file_dir $isoform result";
print "\n$merge_cmd\n";
system($merge_cmd);


$default_score = lc($default_score);
my $prefix = $work_dir."/".$exp_name."Ribo";
my $command_comp = "Rscript $script_dir"."/FDR.R ".$fdr." ".$cutoff." ".$default_score." ".$isoform." ".$prefix;
print "\n$command_comp\n";
system($command_comp);

print "\n";
$msg = "\tAnalysis completed and results written to file.";
print_localtime($msg);
timer($startRun);





###------------- SUBROUTINES -------------###

sub ZINB_parallel {
	my $pm = Parallel::ForkManager->new($cores);
	opendir(DIRECTORY, $file_dir) or die $!;
	while (my $file = readdir(DIRECTORY)) {
		#next if ($file =~ /_pos/);
		my $pid = $pm->start and next;	# parallel process access to split directory
		unless ($file=~/^\./){
			my $file_in = path($file, $file_dir);
			system("Rscript $script_dir"."/RiboZINB.R ".$tmp_dir." ".$file_in." ".$file_in."_zinb ".$threshold." ".$no_of_samples." ".$sigma." ".$default_score);
		}
		$pm->finish; 
	}
	close DIRECTORY;
	$pm->wait_all_children;
}


sub split_genes_file {

	my $transcripts = $_[0];

	my $number_gene = scalar(keys %$transcripts);
	my $genes_per_file = ceil($number_gene/$cores);	# number of genes in each file
	my $fcount = 1;		# Count number of files
	my $countgenes = 1;		# Exit count

	my $file = $file_dir."/split_".$fcount;
	open F, ">".$file  or die $!;
	print F "gene\tgene_name\ttranscript\tstrand\tCCDS\tCCDS_ID\tbiotype\ttsl\tlength_tr\trpkm_tr\tlength_cds\tcoverage_cds\trpkm_cds\treads\n";

	foreach my $gene (keys %$transcripts) {
		if ($countgenes % $genes_per_file == 0) {
			close F;
			$fcount++;
			my $file = $file_dir."/split_".$fcount;

			open F, ">".$file  or die $!;
			print F "gene\tgene_name\ttranscript\tstrand\tCCDS\tCCDS_ID\tbiotype\ttsl\tlength_tr\trpkm_tr\tlength_cds\tcoverage_cds\trpkm_cds\treads\n";
			$countgenes = 1;
		}
	
		foreach my $tr (sort keys %{$transcripts->{$gene}}) {
			my $gene_name = $transcripts->{$gene}->{$tr}->{gene_name};
			my $strand = $transcripts->{$gene}->{$tr}->{strand};
			my $CCDS = $transcripts->{$gene}->{$tr}->{CCDS};
			my $CCDS_ID = $transcripts->{$gene}->{$tr}->{ccds_id};
			my $tr_biotype = $transcripts->{$gene}->{$tr}->{biotype};
			my $rpkm_tr = $transcripts->{$gene}->{$tr}->{rpkm};
			my $coverage_cds = $transcripts->{$gene}->{$tr}->{coverage_cds};
			my $rpkm_cds = $transcripts->{$gene}->{$tr}->{rpkm_cds};
			my $reads_tr = $transcripts->{$gene}->{$tr}->{reads};
			my $length_cds = $transcripts->{$gene}->{$tr}->{length_cds};
			my $length_tr = $transcripts->{$gene}->{$tr}->{len};
 			my $tsl = $transcripts->{$gene}->{$tr}->{tsl};

            unless ($CCDS_ID) {$CCDS_ID = ""}
			print F "$gene\t$gene_name\t$tr\t$strand\t$CCDS\t$CCDS_ID\t$tr_biotype\t$tsl\t$length_tr\t$rpkm_tr\t$length_cds\t$coverage_cds\t$rpkm_cds\t$reads_tr\n";
		}

		$countgenes++;	
	}
}

sub transcript_info {

	my $transcripts = $_[0];
	my $filename 	= $_[1];

	open F, ">".$filename or die $!;
	print F "gene\tgene_name\ttranscript\tstrand\tCCDS\tCCDS_ID\tbiotype\ttsl\tlength_tr\trpkm_tr\tlength_cds\tcoverage_cds\trpkm_cds\treads\n";

	foreach my $gene (sort keys %$transcripts) {
		foreach my $tr (sort keys %{$transcripts->{$gene}}) {
			my $gene_name = $transcripts->{$gene}->{$tr}->{gene_name};
			my $CCDS = $transcripts->{$gene}->{$tr}->{CCDS};
			my $CCDS_ID = $transcripts->{$gene}->{$tr}->{ccds_id};
			my $strand = $transcripts->{$gene}->{$tr}->{strand};
			my $tr_biotype = $transcripts->{$gene}->{$tr}->{biotype};
			my $rpkm_tr = $transcripts->{$gene}->{$tr}->{rpkm};
			my $coverage_cds = $transcripts->{$gene}->{$tr}->{coverage_cds};
			my $rpkm_cds = $transcripts->{$gene}->{$tr}->{rpkm_cds};
			my $reads_tr = $transcripts->{$gene}->{$tr}->{reads};
			my $length_cds = $transcripts->{$gene}->{$tr}->{length_cds};
			my $length_tr = $transcripts->{$gene}->{$tr}->{len};
 			my $tsl = $transcripts->{$gene}->{$tr}->{tsl};

            unless ($CCDS_ID) {$CCDS_ID = ""}
			print F "$gene\t$gene_name\t$tr\t$strand\t$CCDS\t$CCDS_ID\t$tr_biotype\t$tsl\t$length_tr\t$rpkm_tr\t$length_cds\t$coverage_cds\t$rpkm_cds\t$reads_tr\n";
		}
	}
	close F;
}


sub gene_info {

	my $genes  = $_[0];
	my $filename = $_[1];

	open F, ">".$filename or die $!;
	print F "gene\tgene_name\tchromosome\tstrand\tstart\tstop\tbiotype\treads\trpkm\tRPF_pos\tno_of_transcripts\tproportion_of_zero\n";
	foreach my $gene (sort keys %$genes) {
		my $gene_name	= $genes->{$gene}->{gene_name};
		my $chr 	= $genes->{$gene}->{region};
		my $strand 	= $genes->{$gene}->{strand};
		my $biotype = $genes->{$gene}->{biotype};
		my $reads	= $genes->{$gene}->{reads};
		my $rpkm	= $genes->{$gene}->{rpkm};
		my $RPF_pos = $genes->{$gene}->{RPF_pos};
		my $tr_count = $genes->{$gene}->{tr_count};
		my $prop_of_zero = $genes->{$gene}->{prop_of_zero};
		my $start = $genes->{$gene}->{start};
		my $stop = $genes->{$gene}->{stop}; 
		print F "$gene\t$gene_name\t$chr\t$strand\t$start\t$stop\t$biotype\t$reads\t$rpkm\t$RPF_pos\t$tr_count\t$prop_of_zero\n";
	}
	close F;
}


sub Isoform_data {

	# statistical model requires at least 50 data points for accurate prediction
	# Consider only genes with a minimum of MINREADS count positions pad positions that don't exist for each transcript by zero
	my $genes = $_[0];
	my $transcripts = $_[1];

	my $transcripts_valid = {};
	my $valid_genes = {};

	foreach my $gene (sort keys %$transcripts) {
		my $chr = $genes->{$gene}->{region};
		my $strand 	= $genes->{$gene}->{strand};
		my $biotype = $genes->{$gene}->{biotype};

		unless ($biotype) {print "$gene\n"; exit}
		if ($biotype eq 'protein_coding') {
			my $gene_positions = {};	# hash to store gene positional information
			foreach my $exon_id (keys %{$genes->{$gene}->{exons}}){
				my @exon = split '-', $genes->{$gene}->{exons}->{$exon_id};
				for(my $p = $exon[0]; $p <=$exon[1]; $p++) {
					if (exists $reads_count_table->{$chr}->{$strand}->{$p}) {
						$gene_positions->{$p} = $reads_count_table->{$chr}->{$strand}->{$p};
					} else {
						$gene_positions->{$p} =  0;
					}
				}
			}

			my $RPF_pos = 0;
			my $gene_reads = 0;
			foreach my $p (keys %$gene_positions) {
				if ($gene_positions->{$p} > 0) {$RPF_pos++}
				$gene_reads += $gene_positions->{$p};
			}

			# collect information about the genes for downstream analysis
			my @positions = ($strand eq '+') ? sort {$a <=> $b} keys %$gene_positions: sort {$b <=> $a} keys %$gene_positions;
			my $gene_coverage = $RPF_pos/scalar(@positions);
            my $gene_rpkm = ($gene_reads*1000000000)/(scalar(@positions)*$mapped_total);

			$genes->{$gene}->{reads} = $gene_reads;
			$genes->{$gene}->{prop_of_zero} = 1 - $gene_coverage;
			$genes->{$gene}->{RPF_pos} = $RPF_pos;
			$genes->{$gene}->{rpkm} = $gene_rpkm;

			# to ensure the zeros are inflated in all cases
			if ($gene_reads >= $MINCOUNT) {
				my $count_tr = 0;	# count number of valid trnscript per gene 
				my $tr_hash = {};	# hash to store each transcript genomic positions
				my @trans = ();	 	# array to keep track of all transcript of the currest gene

				foreach my $tr (sort keys %{$transcripts->{$gene}}) {

					my $tr_biotype = $transcripts->{$gene}->{$tr}->{biotype};

					if ($tr_biotype eq 'protein_coding' or $tr_biotype eq 'Polymorphic pseudogene') {
						my @exons = @{$transcripts->{$gene}->{$tr}->{exons}};
						my $start_cds = $transcripts->{$gene}->{$tr}->{start_cds};
						my $stop_cds = $transcripts->{$gene}->{$tr}->{stop_cds};
						unless ($start_cds) {$start_cds = $transcripts->{$gene}->{$tr}->{start}}
						unless ($stop_cds) {$stop_cds = $transcripts->{$gene}->{$tr}->{stop}}

						my ($tr_positions,$tr_reads,$rpkm,$coverage_cds,$rpkm_cds,$length_cds) = transcript_positions(\@exons,$chr,$strand,$start_cds,$stop_cds);
						$transcripts->{$gene}->{$tr}->{rpkm} = $rpkm;
						$transcripts->{$gene}->{$tr}->{coverage_cds} = $coverage_cds;
						$transcripts->{$gene}->{$tr}->{rpkm_cds} = $rpkm_cds;
						$transcripts->{$gene}->{$tr}->{reads} = $tr_reads;
						$transcripts->{$gene}->{$tr}->{length_cds} = $length_cds;
						$transcripts->{$gene}->{$tr}->{len} = scalar(keys %{$tr_positions});

						# keep track of all positions in tr_hash
						foreach my $p (keys %$tr_positions) {
							$tr_hash->{$tr}->{$p} = $tr_positions->{$p};
						}
						$count_tr++;
						push @trans, $tr;	# add transcript into the array
						$transcripts_valid->{$gene}->{$tr} = $transcripts->{$gene}->{$tr};
					}
				}

				# allow only gene with at least one transcript
                my $gene_ave_read = ceil($gene_reads/scalar(@positions));
				my $filename = $tmp_dir."/".$gene.".txt";
				open F, ">".$filename or die $!;
				print F "$gene\t",join("\t", @trans);
				for (my $i = 0; $i <scalar(@positions); $i++) {

					my $p = $positions[$i];
					print F "\n$gene_positions->{$p}";

					for (my $t=0; $t<scalar(@trans); $t++) {
					    my $rpf_count = (exists $tr_hash->{$trans[$t]}->{$p}) ? $tr_hash->{$trans[$t]}->{$p}: 0;
						print F "\t$rpf_count";
					}
				}
				close F;

				$valid_genes->{$gene} = $genes->{$gene};
				$valid_genes->{$gene}->{tr_count} = $count_tr;
			}
		}
	}

	if (scalar(keys %$valid_genes) == 0) {
		print "No gene found with RPF reads above the minimum of $MINCOUNT reads\n";
		exit;
	}

	return ($valid_genes,$transcripts_valid);

}


sub transcript_positions {

	my $exons 	 = $_[0];
	my $chr 	 = $_[1];
	my $strand 	 = $_[2];
	my $start_cds= $_[3];
	my $stop_cds = $_[4];

	my $tr_positions = {};	# store reads for each nucleotide within CDS
	my $reads = 0;
	my $rpkm_cds = 0;
	my $coverage_cds = 0;
	my $length_cds = 0;

	foreach my $exon (@{$exons}){
		for(my $p = $exon->[0]; $p <= $exon->[1]; $p++) {

			if ($reads_count_table->{$chr}->{$strand}->{$p}) {
				$tr_positions->{$p} = $reads_count_table->{$chr}->{$strand}->{$p};
				$reads += $reads_count_table->{$chr}->{$strand}->{$p};
				if ($start_cds <= $p and $p <= $stop_cds) {
					$rpkm_cds += $reads_count_table->{$chr}->{$strand}->{$p};
					$coverage_cds++;
				}
			} else {
				$tr_positions->{$p} =  0;
			}

			if ($start_cds <= $p and $p <= $stop_cds) {$length_cds++}
		}
	}

    $length_cds = $length_cds - 1;
	my $rpkm = ($reads*1000000000)/(scalar(keys %$tr_positions)*$mapped_total);
	$rpkm_cds = ($rpkm_cds*1000000000)/($length_cds*$mapped_total);
	$coverage_cds = $coverage_cds/$length_cds;

	return $tr_positions,$reads,$rpkm,$coverage_cds,$rpkm_cds,$length_cds;

}


sub annotation_GTF {

	my $file 	= $_[0];

	my $transcripts = {};
	my $genes = {};
	
	open (F, $file) or die  "Error reading file: $file";
	while (<F>) {
		next if (/^#/);
		chomp $_;

		my @line 	= (split '\t', $_);
		my $chr 	= $line[0];
		my $source= $line[1];
		my $feature= $line[2];
		my $start = $line[3];
		my $stop 	= $line[4];
		my $strand= $line[6];
		my $info	= $line[8];

		$chr =~ s/^chr//;
		if ($feature eq 'gene') {

			my ($gene_id) = $info =~ /gene_id."?([^";]+)"?/;
			my ($gene_name) = $info =~ /gene_name."?([^";]+)"?/;
			my ($biotype) = $info =~ /gene_biotype."?([^";]+)"?/;
			unless ($gene_name) {$gene_name = 'na'}

			$genes->{$gene_id}->{gene_name}= $gene_name;
			$genes->{$gene_id}->{strand}= $strand;
		 	$genes->{$gene_id}->{region}= $chr;
		 	$genes->{$gene_id}->{biotype}= $biotype;
		 	$genes->{$gene_id}->{start}= $start;
		 	$genes->{$gene_id}->{stop}= $stop;
		}

		if ($feature eq 'transcript') {
			my ($tr_id)   = $info =~ /transcript_id."?([^";]+)"?/;
			my ($gene_id) = $info =~ /gene_id."?([^";]+)"?/;
			my ($gene_name) = $info =~ /gene_name."?([^";]+)"?/;
			my ($biotype_tr) = $info =~ /transcript_biotype."?([^";]+)"?/;
			my ($ccds_id) = $info =~ /ccds_id."?([^";]+)"?/;
			my ($tsl) = $info =~ /transcript_support_level."?([^";]+)"?/;
			unless ($gene_name) {$gene_name = 'na'}

			#if ($biotype_tr eq 'protein_coding' or $biotype_tr eq 'Polymorphic pseudogene') {
			if ($biotype_tr eq 'protein_coding' and exists $genes->{$gene_id}) {
		 		$transcripts->{$gene_id}->{$tr_id}->{region} = $chr;
		 		$transcripts->{$gene_id}->{$tr_id}->{gene_name} = $gene_name;
		 		$transcripts->{$gene_id}->{$tr_id}->{start} = $start;
		 		$transcripts->{$gene_id}->{$tr_id}->{stop} 	= $stop;
		 		$transcripts->{$gene_id}->{$tr_id}->{strand}= $strand;
	 			$transcripts->{$gene_id}->{$tr_id}->{region}= $chr;
	 			$transcripts->{$gene_id}->{$tr_id}->{biotype}= $biotype_tr;
	 			$transcripts->{$gene_id}->{$tr_id}->{tsl} = (defined $tsl) ? $tsl: 'na';
	 			$transcripts->{$gene_id}->{$tr_id}->{ccds_id} = $ccds_id;
	 			$transcripts->{$gene_id}->{$tr_id}->{CCDS} = ($info =~ /tag "CCDS";/) ? 'Y': 'N';
			}
		}

		if ($feature eq 'exon') {
			my ($tr_id)   = $info =~ /transcript_id."?([^";]+)"?/;
			my ($exon_id)  = $info =~ /exon_id."?([^";]+)"?/;
			my ($gene_id) = $info =~ /gene_id."?([^";]+)"?/;
			my ($biotype_tr) = $info =~ /transcript_biotype."?([^";]+)"?/;
			my @exon = ($start,$stop);

			#if ($biotype_tr eq 'protein_coding' or $biotype_tr eq 'Polymorphic pseudogene') {
			if ($biotype_tr eq 'protein_coding' and exists $genes->{$gene_id}) {
 				push @{$transcripts->{$gene_id}->{$tr_id}->{exons}}, \@exon;
				$genes->{$gene_id}->{exons}->{$exon_id} = $start."-".$stop;
			}
		}

		if ($feature eq 'start_codon') {
			my ($tr_id)   = $info =~ /transcript_id."?([^";]+)"?/;
			my ($gene_id) = $info =~ /gene_id."?([^";]+)"?/;
			my ($biotype_tr) = $info =~ /transcript_biotype."?([^";]+)"?/;

			#if ($biotype_tr eq 'protein_coding' or $biotype_tr eq 'Polymorphic pseudogene') {
			if ($biotype_tr eq 'protein_coding' and exists $genes->{$gene_id}) {
		        if ($strand eq '+') {
					$transcripts->{$gene_id}->{$tr_id}->{start_cds} = $start;
		        } else {
					$transcripts->{$gene_id}->{$tr_id}->{stop_cds} = $stop;
		        }
			}
		}

		if ($feature eq 'stop_codon') {
			my ($tr_id)   = $info =~ /transcript_id."?([^";]+)"?/;
			my ($gene_id) = $info =~ /gene_id."?([^";]+)"?/;
			my ($biotype_tr) = $info =~ /transcript_biotype."?([^";]+)"?/;

			#if ($biotype_tr eq 'protein_coding' or $biotype_tr eq 'Polymorphic pseudogene') {
			if ($biotype_tr eq 'protein_coding' and exists $genes->{$gene_id}) {
		        if ($strand eq '+') {
					$transcripts->{$gene_id}->{$tr_id}->{stop_cds} = $start;
		        } else {
					$transcripts->{$gene_id}->{$tr_id}->{start_cds} = $stop;
		        }
			}
		}
	}
	close F;

	return $genes,$transcripts;
}


sub path {
	my ($file, $dir) = @_;
	return (File::Spec->catfile( $dir, $file));
}


sub get_read_table {

	my $file = $_[0];
	my $read = {};		# hash to store all ribo-seq inofrmation from SQLite DB

	open(F, "$file") or die "Cannot open file $file\n";
	while (<F>) {
		chomp $_;
		# chr,strand,start,count
		next if (/^chr/);
		my @line 	= split ',|\t', $_;
		my $chr 	= $line[0];
		my $strand = ($line[1] == 1) ? '+': '-';
		#my $strand = $line[1];
		my $start = $line[2];
		my $count = $line[3];

		$count =~ s/\n|\r//g;
		$chr =~ s/^chr//;
		$read->{$chr}->{$strand}->{$start} = $count;
	}
	
	if (scalar($read) == 0) {
		print "No valid information from in RPF counts per genomic position file. Please ensure the file is in the accepted format [chr,strand,position,RPF_count]\n";
		exit;
	}

	return $read;
}

sub uninitialized_param {
	my ($v) = @_;
	not ( defined($v) and length $v );
}


sub print_localtime {
	my $message = $_[0];
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	printf("%02d:%02d:%02d", $hour, $min, $sec);
	print "$message\n";

}


sub check_work_dir {

	my $work_dir = $_[0];

	my $tmp_dir = $work_dir."/tmp";
	my $file_dir = $work_dir."/gene_split";

	if ($work_dir) {
		if (!-d $work_dir) {
			system("mkdir -p $work_dir");
			system("mkdir -p $tmp_dir");	# recreate tmp directory
			system("mkdir -p $file_dir");	# recreate tmp directory
		} else {
			if (uc($force) eq 'Y') {
				system("rm -rf $work_dir");	# delete content of exiting tmp directory
				system("mkdir -p $work_dir");	# delete content of exiting tmp directory
				system("mkdir -p $tmp_dir");	# recreate tmp directory
				system("mkdir -p $file_dir");	# recreate tmp directory
			} else {
                print "WARNING: The output directory '$work_dir' already exists. \nExisting files will be overwritten.\n";
				if (-d $tmp_dir) {
                    print "WARNING: The output directory '$tmp_dir' already exists. \nExisting files will be overwritten.\n";
				} else {
					system("mkdir -p $tmp_dir");	# recreate tmp directory
				}

				if (-d $file_dir) {
                    print "WARNING: The output directory '$file_dir' already exists. \nExisting files will be overwritten.\n";
				} else {
					system("mkdir -p $file_dir");	# recreate tmp directory
				}
			}
		}
	} else {
		# get date and time to create tmp folder
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
		$work_dir = "RiboZinB_iso".$mday.$mon.$year;

		if (!-d $work_dir) {
			$tmp_dir = $work_dir."/tmp";
			$file_dir = $work_dir."/gene_split";
			system("mkdir -p ".$work_dir);
			system("mkdir -p ".$tmp_dir);	# recreate tmp directory
			system("mkdir -p ".$file_dir);	# recreate tmp directory
		} else {
			$work_dir = "Isoform_prediction".$mday.$mon.$year."_2";
			$tmp_dir = $work_dir."/tmp";
			$file_dir = $work_dir."/gene_split";
			system("mkdir -p ".$work_dir);
			system("mkdir -p ".$tmp_dir);	# recreate tmp directory
			system("mkdir -p ".$file_dir);	# recreate tmp directory
		}
	}

	unless ($file_dir) {print "Directories $file_dir cannot be created ensure you have the appropriate permission\n";exit}
	unless ($tmp_dir) {print "Directories $tmp_dir cannot be created ensure you have the appropriate permission\n";exit}

	return $work_dir, $file_dir, $tmp_dir;
}


sub _isNumeric {
    # check if variable is numeric
    my ($value) = @_; 

    no warnings;

    return 1 if ($value + 0) eq $value;
    return 0;
}

sub timer {
	my $startRun = shift;
	my $endRun 	= time();
	my $runTime = $endRun - $startRun;
	printf("\nTotal running time: %02d:%02d:%02d\n\n", int($runTime / 3600), int(($runTime  % 3600) / 60), int($runTime % 60));
}


