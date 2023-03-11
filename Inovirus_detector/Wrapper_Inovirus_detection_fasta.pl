#!/usr/bin/env perl
use strict;
use autodie;
use File::Basename;
use Bio::SeqIO;
use Getopt::Long;
use Cwd;
my $h='';
my $cmd='';
my $out='';
my $original_fna_file='';
my $path_pfam='';
my $n_cpu=2;
my $path_signalp='';
my $path_signalp5='';
my $path_tmhmm='';
my $db_dir="Inovirus_db/";
my $tag_second_only=0;
my $min_score=0.83;
GetOptions ('help' => \$h, 'h' => \$h, 'f=s'=>\$original_fna_file, 'p=s'=>\$path_pfam, 'd=s'=>\$db_dir, 't=s'=>\$n_cpu, 'sp5=s'=>\$path_signalp5 , 'sp=s'=>\$path_signalp, 'th=s'=>\$path_tmhmm, 's'=>\$tag_second_only, 'min_score=s'=>\$min_score);
if ($h==1 || $original_fna_file eq "" || $path_pfam eq "" || $path_tmhmm eq "" || ($path_signalp eq "" && $path_signalp5 eq "")){ # If asked for help or did not set up any argument
	print "# Script to predict putative inovirus sequences from a gb file
#### Arguments :
# -f : fna file of the contigs to be analyzed
# -p : path to the pfam database hmm collection, i.e. Pfam-A.hmm (ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz)
# -th : path to tmhmm (http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm)
# -sp : path to signal_p <v4 (http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp)
# or
# -sp5 : path to signal_p v5 folder (www.cbs.dtu.dk/services/SignalP/portable.php)
#### Optional arguments
# -d : path to the Inovirus_db folder (default: Inovirus_db/)
# -t : number of threads used for hmmsearch and blast (default: 2)
# -s : run only the second step
# -min_score: change the minimum score cutoff (default = 0.83)
#### Requirements:
# Bio::SeqIO
# Hmmsearch
# blastp
# prodigal
# aragorn
# Predict_inovirus_coat_proteins.pl
####
";
	die "\n";
}

if ($tag_second_only==1){
	print "Skipping step 1\n";
}
else{
	&run_cmd("./Identify_candidate_fragments_from_fna.pl -f $original_fna_file -p $path_pfam -th $path_tmhmm -sp $path_signalp"); ## Note - need to adjust for Sp5
}

my $tag=0;
if ($original_fna_file=~/(.*)\/([^\/]+)$/){
	my $dir=$1;
	my $suffix=$2;
	if ($suffix=~/(.*)\.[^\.]+$/){$suffix=$1;}
	my $out_summary=$dir."/".$suffix."_prediction_summary.tsv";
	foreach my $annot_file (<$dir/*annot.gff>){
		my $fna_file=$annot_file;
		$fna_file=~s/_annot.gff$/_nucl.fna/;
		&run_cmd("./Get_inovirus_prediction_score_from_gff_fragments.pl -i $annot_file -f $fna_file -th $min_score");
	}
	## Prepare the summary
	$tag=0;
	open my $s1,">",$out_summary;
	print $s1 "Fragment_id\tScore\tDetection type\tStart\tStop\tLength\tDR identity %\tDR length\n";
	foreach my $annot_file (<$dir/*annot.gff>){
		$annot_file=~/$dir\/(.*)_annot.gff/;
		my $frag_id=$1;
		my $test_1=$dir."/".$frag_id."_annot_inovirus-predictions.csv";
		my $test_2=$dir."/".$frag_id."_annot_inovirus-predictions-refined.csv";
		if (-e $test_1){
			if (-e $test_2){
				$tag=0;
				open my $csv,"<",$test_2;
				while(<$csv>){
					chomp($_);
					if ($_=~/^#/){next;}
					if ($_=~/^>(.*)/){
						my $line=$1;
						my @tab=split(",",$line);
						my @t=split("_",$tab[1]);
						my @t1=split("-",$t[$#t]);
						if ($tab[5] eq ""){$tab[5]="NA";}
						if ($tab[6] eq ""){$tab[6]="NA";}
						print $s1 $tab[1]."\t".$tab[4]."\t".$tab[0]."\t".$t1[0]."\t".$t1[1]."\t".$tab[2]."\t".$tab[5]."\t".$tab[6]."\n";
						$tag=1;
					}
				}
				close $csv;
				# if ($tag==0){
				# 	print $s1 $frag_id."\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
				# }
			}
			else{
				$tag=0;
				open my $csv,"<",$test_1;
				while(<$csv>){
					chomp($_);
					if ($_=~/^#/){next;}
					if ($_=~/^>(.*)/){
						my $line=$1;
						my @tab=split(",",$line);
						my @t=split("_",$tab[1]);
						my @t1=split("-",$t[$#t]);
						if ($tab[5] eq ""){$tab[5]="NA";}
						if ($tab[6] eq ""){$tab[6]="NA";}
						print $s1 $tab[1]."\t".$tab[4]."\t".$tab[0]."\t".$t1[0]."\t".$t1[1]."\t".$tab[2]."\t".$tab[5]."\t".$tab[6]."\n";
						$tag=1;
					}
				}
				close $csv;
				# if ($tag==0){
				# 	print $s1 $frag_id."\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
				# }
			}
		}
	}
	close $s1;
}
else{
	print "I do not understand $original_fna_file format, so I won't run the second step or the summary, sorry\n";
}






sub run_cmd{
	my $cmd=$_[0];
	if ($_[1] ne "veryquiet"){print "$cmd\n";}
	if ($_[1] ne "dummy"){
		my $out=`$cmd`;
		if ($_[1] ne "quiet" && $_[1] ne "veryquiet"){
			if ($_[1] eq "stderr"){print STDERR "$out\n";}
			else{print "$out\n";}
		}
		return($out);
	}
	else{
		print " ### dummy run\n";
	}
}
