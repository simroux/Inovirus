#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h='';
my $cmd='';
my $out='';
my $fa_file='';
my $gff_file='';
my $th=0.83; ## Suggested (empirical) threshold
my $db_dir="Inovirus_db/";
my $path_run_predict="run_predict.R";
my $path_r_model="InovirusModel.RData";
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$gff_file, 'f=s'=>\$fa_file, 'd=s'=>\$db_dir, 'prt=s'=>\$path_run_predict, 'prm=s'=>\$path_r_model, 'th=s'=>\$th);
if ($h==1 || $gff_file eq "" || $fa_file eq ""){ # If asked for help or did not set up any argument
	print "# Script to get the an Inovirus prediction from contigs / fragments where a putative morphogenesis gene was identified
#### Recommendation is to take a fragment of ~ +/- 30kb around the putative morphogenesis gene
#### Arguments :
# -i : gff of the fragments to be analyzed
# -f : nucleotide fasta file of the fragments to be analyzed
#### Optional arguments
# -d : path to the Inovirus_db folder (default: Inovirus_db/)
# -prt : path to the R script 'run_predict.R' (default: run_predict.R)
# -prm : path to the R model 'InovirusModel.RData' (default: InovirusModel.RData)
# -th: custom cutoff on model score (default: 0.83, should be between 0 and 1)
#### Requirements:
# Blast+
# R
# InovirusModel.RData (which includes the custom sequence classifier)
#### GFF format:
# Each fragment should be preceeded by a line starting with \"## Fragment_id\", The Fragment_id must be the gene id of the putative morphogenesis gene
# PFAM annotation should be indicated as \"pfam=PFAM_domain_name Score\" in the 9th field
# Inoviridae PC annotation with score should be indicated as \"inoPC=PC_XXX Score\" in the 9th field as well
# Prediction of coat protein should be indicated as \"Coat_pred=yes/no\" in the 9th field as well
# Note: the fragment is expected to be a microbial genome fragment, and start with a CDS.
# The starting coordinate of the first CDS will thus be used to \"shift\" the coordinate when dealing with the fragment sequence (i.e. transform the genes coordinates from the global genome reference to the fragment reference)
# if it's not the case, you can indicate in the header line after the fragment id (tab separated) a number of bp to use in this transformation (can be 0 if the fragment is an entire contig)
#### Fasta format:
# Each fragment should have its nucleotide sequence in the fasta file, identified by the same Fragment_id as in the gff file (i.e. gene id of the putative morphogenesis gene
####
";
	die "\n";
}

if (!($db_dir=~/\/$/)){$db_dir.="/";}
if (!(-d $db_dir)){die("pblm, we did not find the directory Inovirus_db -> $db_dir ?\n");}

my $expected_pfam_list=$db_dir."List_known_PFAM.tsv";
my $nonduf_list=$db_dir."List_DUF_affiliation.txt";


my $test_path=&run_cmd('which Rscript','quiet');
if ($test_path eq ""){die("pblm, we did not find a path for Rscript, was the module/conda environment loaded ? ($test_path)\n");}


$gff_file=~/(.*)\.[^\/]+/;
my $basename=$1;
my $final_out=$basename."_inovirus-predictions.csv";
my $final_out_2=$basename."_inovirus-predictions-refined.csv";
my $wdir="./";
if ($basename=~/(.*\/)[^\/]+/){$wdir=$1;}
print "Wdir: ".$wdir."\n";


my %expected;
open my $tsv,"<",$expected_pfam_list;
while(<$tsv>){
	chomp($_);
	if ($_=~/^#/ || $_ eq ""){next;}
	my @tab=split("\t",$_);
	$expected{$tab[0]}=1;
}
close $tsv;

my %check_nonduf;
open my $tsv,"<",$nonduf_list;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	$check_nonduf{$tab[0]}=1;
}
close $tsv;


my %check_frag;
my $c_c="";
open my $fa,"<",$fa_file;
while(<$fa>){
	chomp($_);
	if ($_=~/^>(\S+)/){
		$c_c=$1;
		print "Taking sequence and length for fragment $c_c\n";
	}
	else{
		$check_frag{$c_c}{"contig_length"}+=length($_);
		$check_frag{$c_c}{"seq"}.=$_;
	}
}
close $fa;




## Load the annotation of each gene
my %store_gene;
my %shift;
open my $csv,"<",$gff_file;
while(<$csv>){
	chomp($_);
	if ($_=~/^## (.*)/){
		my @t=split("\t",$1);
		$c_c=$t[0];
		if ($t[1] > -1 ){
			$shift{$c_c}=$t[1];
		}
		print "Reading genes for fragment $c_c\n";
		next;
	}
	my @tab=split("\t",$_);
	my $gene_id="";
	if ($tab[8]=~/ID=([^;]+)/){$gene_id=$1;}
	if ($gene_id eq ""){next;}
	if (!defined($shift{$c_c})){$shift{$c_c}=$tab[3];}
	$store_gene{$c_c}{$gene_id}{"start"}=$tab[3];
	$store_gene{$c_c}{$gene_id}{"start_shifted"}=$tab[3]-$shift{$c_c};
	$store_gene{$c_c}{$gene_id}{"stop"}=$tab[4];
	$store_gene{$c_c}{$gene_id}{"stop_shifted"}=$tab[4]-$shift{$c_c};
	$store_gene{$c_c}{$gene_id}{"strand"}=$tab[6];
	$store_gene{$c_c}{$gene_id}{"type"}=$tab[2];
	$store_gene{$c_c}{$gene_id}{"line"}=$tab[0].",".$tab[3].",".$tab[4].",".$store_gene{$c_c}{$gene_id}{"start_shifted"}.",".$store_gene{$c_c}{$gene_id}{"stop_shifted"}.",".$tab[6].",".$tab[2].",".$tab[8];
	## Load Pfam affiliation
	if ($tab[8]=~/pfam=([^;]+)/){
		if ($1 ne "NA"){
			my @t=split(" ",$1);
			if ($t[1]>=30){
				$store_gene{$c_c}{$gene_id}{"pfam"}=$t[0];
				print "$gene_id -> $t[0] - PFAM\n";
# 				<STDIN>;
			}
		}
	}
	elsif($tab[2] eq "CDS"){
		print "?! No Pfam affiliation for $gene_id ?\n";
		<STDIN>;
	}
	## Load inovirus PC affiliation
	if ($tab[8]=~/inoPC=([^;]+)/){
		if ($1 ne "NA"){
			my @t=split(" ",$1);
			if ($t[1]>=30){
				$store_gene{$c_c}{$gene_id}{"pc"}=$t[0];
				print "$gene_id -> $t[0] - InoPC\n";
			}
		}
	}
	elsif($tab[2] eq "CDS"){
		print "?! No Inovirus affiliation for $gene_id ?\n";
		<STDIN>;
	}
	## Load Coat protein prediction
	if ($tab[8]=~/coat_pred=([^;]+)/){
		if ($1 ne "NA"){
			$store_gene{$c_c}{$gene_id}{"pred"}="coat";
			print "$gene_id -> $1 - CoatPred\n";
		}
	}
	else{
		print "?! No Coat prediction for $gene_id ?\n";
		<STDIN>;
	}
}
close $csv;



my $tag_detection=0;
my %info_frag;
open my $s_out,">",$final_out;
print $s_out "## >Fragment type,Fragment Id_start-stop,final_fragment_length,initial_fragment_length,Prediction_confidence_score\n";
print $s_out "## Contig,Start,Stop,Start_transformed,Stop_transformed,Strand,Type,Affiliations\n";
foreach my $frag (sort keys %check_frag){
	print "Working on fragment $frag\n";
	# Loading all genes from this genome
	my @tab_genes=sort {$store_gene{$frag}{$a}{"start"} <=> $store_gene{$frag}{$b}{"start"}} keys %{$store_gene{$frag}};
	my $tmp_input=$wdir."Tmp_input.csv";
	my %store_result;
	my $n_ini=-1;
	my $nb_genes=0;
	# Check which gene is the putative morphogenesis one
	for (my $i=0;$i<=$#tab_genes;$i++){
		$nb_genes++;
		if (defined($check_frag{$tab_genes[$i]})){$n_ini=$i;}
# 		print "Checking $tab_genes[$i]\n";
	}
	my $start_i=0;
	my $end_i=$#tab_genes;
	open my $s1,">",$tmp_input;
	print $s1 "Id,Type,Length,N_genes,N_morph,N_pfam_coat,N_pred_coat,N_ino_pcs,N_unpfam,R_unpfam,L10,L50\n";
	for (my $n_start=$start_i;$n_start<=$n_ini;$n_start++){
		for (my $n_end=$n_ini;$n_end<=$end_i;$n_end++){
			my $start=$store_gene{$frag}{$tab_genes[$n_start]}{"start"};
			my $end=$store_gene{$frag}{$tab_genes[$n_end]}{"stop"};
			my $id_frag="Temp_".$frag."_".$start."-".$end;
			print "Checking fragment $n_start to $n_end => $start to $end => $id_frag \n";
			my $length_bp=$end-$start;
			my %count=("morph"=>0,"pfam_coat"=>0,"pred"=>0,"pc"=>0,"unpfam"=>0); # initialize the counts to 0
			my @tab_size=();
			my $total_genes=0;
			foreach my $gene (@tab_genes){
				if ($store_gene{$frag}{$gene}{"start"}>=$start && $store_gene{$frag}{$gene}{"start"}<=$end && $store_gene{$frag}{$gene}{"stop"}>=$start && $store_gene{$frag}{$gene}{"stop"}<=$end){
					print "$gene ".$store_gene{$frag}{$gene}{"start"}." - ".$store_gene{$frag}{$gene}{"stop"}." is in the fragment ".$store_gene{$frag}{$gene}{"pfam"}." // ".$store_gene{$frag}{$gene}{"pc"}." // ".$store_gene{$frag}{$gene}{"pred"}."\n";
# 					<STDIN>;
					$total_genes++;
					my $size=($store_gene{$frag}{$gene}{"stop"}-$store_gene{$frag}{$gene}{"start"}+1)/3;
					push(@tab_size,$size);
					if ($store_gene{$frag}{$gene}{"pfam"} ne ""){
						if ($store_gene{$frag}{$gene}{"pfam"}=~/coat/i || $store_gene{$frag}{$gene}{"pfam"}=~/capsid/i){
							$count{"pfam_coat"}++;
							print "Found ".$store_gene{$frag}{$gene}{"pfam"}." as a capsid PFAM domain\n";
# 							<STDIN>;
						}
						if ($expected{$store_gene{$frag}{$gene}{"pfam"}}==1){
							print $store_gene{$frag}{$gene}{"pfam"}." is expected\n";
# 							<STDIN>;
						}
						else{
# 							print $store_gene{$frag}{$gene}{"pfam"}." is unexpected\n";
							if ($store_gene{$frag}{$gene}{"pfam"}=~/DUF/i || $store_gene{$frag}{$gene}{"pfam"}=~/HTH/i || $store_gene{$frag}{$gene}{"pfam"}=~/DNA/i || $store_gene{$frag}{$gene}{"pfam"}=~/repeat/i || $store_gene{$frag}{$gene}{"pfam"}=~/toxin/i || $store_gene{$frag}{$gene}{"pfam"}=~/regul/){
								if ($check_nonduf{$store_gene{$frag}{$gene}{"pfam"}}==1){
									$count{"unpfam"}++;
								}
								else{
									print $store_gene{$frag}{$gene}{"pfam"}." is not unexpected\n";
								}
							}
							else{
								$count{"unpfam"}++;
							}
						}
					}
# 					if ($store_gene{$frag}{$gene}{"type"} ne "CDS"){
# 						$count{"unpfam"}++;
# 					}
					if ($store_gene{$frag}{$gene}{"pc"} ne "" || $store_gene{$frag}{$gene}{"pfam"}=~/Zot/){
						$count{"pc"}++;
						if ($store_gene{$frag}{$gene}{"pc"}=~/PC_006/ || $store_gene{$frag}{$gene}{"pc"}=~/PC_055/ || $store_gene{$frag}{$gene}{"pc"}=~/PC_036/ || $store_gene{$frag}{$gene}{"pfam"}=~/Zot/){
							$count{"morph"}++;
							print "found $gene -> ".$store_gene{$frag}{$gene}{"pc"}." as a putative morphogenesis\n";
# 							<STDIN>;
						}
					}
					if ($store_gene{$frag}{$gene}{"pred"} ne ""){
						$count{"pred"}++;
						print "found $gene as a putative coat protein\n";
					}
# 					<STDIN>;
				}
			}
			@tab_size=sort {$a <=> $b} @tab_size;
			# Unpredicted PFAMs list, N Zot, N Coat_PFAM, N Coat_predicted, N PCs, N non-PFAMs, % non-PFAMs, Gene 1p, Gene 2p, Gene 3p, Gene 4p, Median
			my $ratio=0;
			if ($total_genes>0){$ratio=$count{"unpfam"}/$total_genes;}
			else{$total_genes=0;}

			my $line=$id_frag.",Unsure,".$length_bp.",".$total_genes.",".$count{"morph"}.",".$count{"pfam_coat"}.",".$count{"pred"}.",".$count{"pc"}.",".$count{"unpfam"}.",".$ratio.",".$tab_size[int($#tab_size/10)].",".$tab_size[int($#tab_size/2)];
			print $s1 $line."\n";
			$store_result{$id_frag}{"n_genes"}=$total_genes;
		}
	}
	close $s1;
	# Now do the prediction
	my $tmp_output=$wdir."Tmp_output.csv";
	&run_cmd("Rscript $path_run_predict $tmp_input $tmp_output $path_r_model");
	# Take the largest window with score > 0.9 (if any)
	open my $csv,"<",$tmp_output;
	while(<$csv>){
		chomp($_);
		$_=~s/\"//g;
		my @tab=split(",",$_);
		if ($tab[1]=~/.*_(\d+)-(\d+)/){$store_result{$tab[1]}{"length"}=$2-$1;}
		else{next;}
		$store_result{$tab[1]}{"value"}=$tab[3];
# 		if ($tab[3]>0.7){
# 			print $tab[1]." -> ".$tab[3]."\n";
# 		}
	}
	close $csv;
	my %final;
	$final{"max_noth"}=0;
	foreach my $frag (sort {$store_result{$b}{"n_genes"} <=> $store_result{$a}{"n_genes"}} keys %store_result){
		if ($final{"max_noth"}<$store_result{$frag}{"value"}){$final{"max_noth"}=$store_result{$frag}{"value"}}
		if ($store_result{$frag}{"value"}>=$th){
			if (!defined($final{"max"}) || $store_result{$frag}{"value"}>$final{"max"}{"value"}){
				$final{"max"}{"value"}=$store_result{$frag}{"value"};
				$final{"max"}{"n_genes"}=$store_result{$frag}{"n_genes"};
				$final{"max"}{"frag"}=$frag;
			}
			if (!defined($final{"n_genes"}) || $store_result{$frag}{"n_genes"}>$final{"max"}{"n_genes"}){
				$final{"n_genes"}{"value"}=$store_result{$frag}{"value"};
				$final{"n_genes"}{"n_genes"}=$store_result{$frag}{"n_genes"};
				$final{"n_genes"}{"frag"}=$frag;
			}
		}
	}
	print "Max proba: ".$final{"max_noth"}."\n";
	print "Final max proba above th of $th: ".$final{"max"}{"frag"}.",".$final{"max"}{"n_genes"}.",".$final{"max"}{"value"}."\n";
	print "Final n_genes: ".$final{"n_genes"}{"frag"}.",".$final{"n_genes"}{"n_genes"}.",".$final{"n_genes"}{"value"}."\n";

	if ($final{"max"}{"frag"} eq ""){
		print "## No inovirus prophage/sequence detected in fragment $frag\n";
	}
	else{
		$tag_detection=1;
		## We like max_proba, even if conservative
		# Once we have our likely fragment, we trim it to remove everything after a non-CDS (e.g. a tRNA or RNA)
		my $n_start_final=-1;
		my $n_stop_final=-1;
		$final{"max"}{"frag"}=~/.*_(\d+)-(\d+)/;
		my $start=$1;
		my $end=$2;
		my $n_gene_frag=0;
		for (my $i=0;$i<=$#tab_genes;$i++){
			my $gene=$tab_genes[$i];
			if ($store_gene{$frag}{$gene}{"start"}>=$start && $store_gene{$frag}{$gene}{"start"}<=$end && $store_gene{$frag}{$gene}{"stop"}>=$start && $store_gene{$frag}{$gene}{"stop"}<=$end){
				$n_gene_frag++;
				if ($store_gene{$frag}{$gene}{"start"}==$start && $n_start_final==-1){$n_start_final=$i}
				if ($store_gene{$frag}{$gene}{"stop"}==$end && $n_stop_final==-1){$n_stop_final=$i; print "$i\n";}
	# 			print "$gene ".$store_gene{$frag}{$gene}{"start"}." - ".$store_gene{$frag}{$gene}{"stop"}." is in the fragment \n";
				if ($store_gene{$frag}{$gene}{"type"} ne "CDS"){
					print "### $gene is not a CDS\n";
					if ($i<$n_ini){
						print "\t so we remove it and all genes before\n";
						$n_start_final=$i+1;
					}
					else{
						print "\t so we remove it and all genes after\n";
						$n_stop_final=$i-1;
					}
				}
			}
		}
	# 	print "$n_start_final to $n_stop_final\n";
		## So we have our final fragment, we print it
		my $type="Prophage";
		my $start=$store_gene{$frag}{$tab_genes[$n_start_final]}{"start"};
		my $end=$store_gene{$frag}{$tab_genes[$n_stop_final]}{"stop"};
		my $frag_length=$end-$start+1;
		if ($frag_length>0.8*$check_frag{$frag}{"contig_length"} || $nb_genes==$n_gene_frag){
			$type="Complete";
			$frag_length=$check_frag{$frag}{"contig_length"};
			$start=0;
			$end=$check_frag{$frag}{"contig_length"};
		}
		print $frag." ".$start."-".$end." ".$frag_length." ".$check_frag{$frag}{"contig_length"}." == ".$nb_genes."-".$n_gene_frag." = ".$type."\n";
		print $s_out ">".$type.",".$frag."_".$start."-".$end.",".$frag_length.",".$check_frag{$frag}{"contig_length"}.",".$final{"max"}{"value"}."\n";
		$info_frag{$frag}{"start"}=$start;
		$info_frag{$frag}{"end"}=$end;
		$info_frag{$frag}{"type"}=$type;
		$info_frag{$frag}{"confidence_score"}=$final{"max"}{"value"};
		$info_frag{$frag}{"data"}=">".$type.",".$frag."_".$start."-".$end.",".$frag_length.",".$check_frag{$frag}{"contig_length"}.",".$final{"max"}{"value"}."\n";
		for (my $i=0;$i<=$#tab_genes;$i++){
			my $gene=$tab_genes[$i];
			if ($store_gene{$frag}{$gene}{"start"}>=$start && $store_gene{$frag}{$gene}{"start"}<=$end && $store_gene{$frag}{$gene}{"stop"}>=$start && $store_gene{$frag}{$gene}{"stop"}<=$end){
				print $s_out $store_gene{$frag}{$gene}{"line"}."\n";
				$info_frag{$frag}{"data"}.=$store_gene{$frag}{$gene}{"line"}."\n";
			}
		}
# 		<STDIN>;
	}
}
close $s_out;

if ($tag_detection==0){
	die("We stop here, no detection\n");
}


print "Now refining these predicted fragments -> looking for att sites for prophages\n";

# Tmp files and global variables
my $db_tmp=$wdir."Db_temp.fasta";
my $query_tmp=$wdir."Query_temp.fasta";
my $out_tmp=$wdir."Blastmp.out";
my $th_len=-1;
my $tag=0;


open my $s_out,">",$final_out_2;
print $s_out "## >Fragment type,Fragment Id_start-stop,final_fragment_length,initial_fragment_length,Prediction_confidence_score,DR identity %,DR length (last 2 columns only if an att site with a Direct Repeat was identified)\n";
print $s_out "## Contig,Start,Stop,Start_transformed,Stop_transformed,Strand,Type,Affiliations\n";
my $k=0;
## Take intergenic around the predictions -> Look for direct repeats outside of these fragments
foreach my $frag (sort keys %info_frag){
	$k++;
	my $start=$info_frag{$frag}{"start"}-$shift{$frag};
	my $end=$info_frag{$frag}{"end"}-$shift{$frag};
	my $seq_c=$check_frag{$frag}{"seq"};
	my %store_att;
	print "Looking at frag $k == $frag => $start, $end\n";
	if ($info_frag{$frag}{"type"} eq "Complete"){
		print "\t This was predicted as a complete contig, we don't look for att site\n";
		# So we don't change anything, and just copy the same information as in the previous file
		print $s_out $info_frag{$frag}{"data"};
		next;
	}
	## Look if tRNA or Integrase
	# Prep the whole fragment: i.e. the predicted prophage +/- 5kb
	my $start_f=$start-5000;
	if ($start_f<0){$start_f=0;}
	my $stop_f=$end+5000;
	if ($stop_f>length($seq_c)){$stop_f=length($seq_c);}
	print "We will look at the fragment from $start_f to $stop_f\n";
	my $seq_frag=substr($seq_c,$start_f,($stop_f-$start_f));
	open my $s1,">",$db_tmp;
	print $s1 ">Fragment\n$seq_frag\n";
	close $s1;
	&run_cmd("makeblastdb -in $db_tmp -title $db_tmp -dbtype nucl");
	# Check that the pI gene is absolutely required
	my %required;
	$required{"start"}=$store_gene{$frag}{$frag}{"start_shifted"};
	$required{"stop"}=$store_gene{$frag}{$frag}{"stop_shifted"};
	print "We'll require the gene betwen $required{start} and $required{stop} (pI-like)\n";
	#
	my $bottom_index=0;
	my $top_index=-1;
	my @tab_genes;
	my $tag_putative_att=0;
	my $tag_canonical=0;
	my $n=0;
	foreach my $gene (sort {$store_gene{$frag}{$a}{"start_shifted"} <=> $store_gene{$frag}{$b}{"start_shifted"}} keys %{$store_gene{$frag}}){
		print "$gene - $store_gene{$frag}{$gene}{type} - $store_gene{$frag}{$gene}{start} - $end - $n\n";
		$store_gene{$frag}{$gene}{"index"}=$n;
		push(@tab_genes,$gene);
		if ($store_gene{$frag}{$gene}{"start_shifted"}==$start){$bottom_index=$n;}
		if ($store_gene{$frag}{$gene}{"stop_shifted"}==$end){$top_index=$n;}
		$n++;
		if ($store_gene{$frag}{$gene}{"type"} eq "tRNA"){$tag_putative_att=1;}
		if ($store_gene{$frag}{$gene}{"pfam"}=~/integrase/i || $store_gene{$frag}{$gene}{"pfam"}=~/^rve/ || $store_gene{$frag}{$gene}{"pfam"}=~/^Transposase_20/){$tag_putative_att=1;}
	}
	print "Top index $top_index -- Bottom index $bottom_index ==== $start == $end\n";
	if ($tag_putative_att==0){
		print "No possibility of a canonical att site here\n";
	}
	else{
		# list tRNAs
		print "Checking possible att sites at tRNAs\n";
		my %putative_att;
		my $mean=($bottom_index+$top_index)/2;
		for (my $i=$bottom_index-5;$i<=$top_index+5;$i++){
			if ($i<0 || $i>$#tab_genes){next;}
			print "$i - ".$store_gene{$frag}{$tab_genes[$i]}{"start"}." - ".$store_gene{$frag}{$tab_genes[$i]}{"stop"}." - ".$store_gene{$frag}{$tab_genes[$i]}{"type"}." - ".$store_gene{$frag}{$tab_genes[$i]}{"pfam"}."\n";
			if (defined($tab_genes[$i])){
				if ($store_gene{$frag}{$tab_genes[$i]}{"type"} eq "tRNA"){
					$putative_att{$i}{"order"}=abs($mean-$i);
				}
			}
		}
		foreach my $index (sort {$putative_att{$a}{"order"} <=> $putative_att{$b}{"order"}} keys %putative_att){
			my $gene=$tab_genes[$index];
			print "Checking gene $gene ($index)\n";
			print "Shift ".$store_gene{$frag}{$gene}{"start_shifted"}." --- ".$store_gene{$frag}{$gene}{"stop_shifted"}."\n";
			my $seq=substr($seq_c,$store_gene{$frag}{$gene}{"start_shifted"},$store_gene{$frag}{$gene}{"stop_shifted"}-$store_gene{$frag}{$gene}{"start_shifted"});
			open my $s1,">",$query_tmp;
			print $s1 ">Query\n$seq\n";
			close $s1;
			&run_cmd("blastn -query $query_tmp -db $db_tmp -out $out_tmp -task blastn-short -num_threads 1 -outfmt 6 -evalue 0.1");
			# We want hits with at least 90% ANI and an att site of 20bp
			my $start_e=$store_gene{$frag}{$gene}{"start_shifted"};
			my $stop_e=$store_gene{$frag}{$gene}{"stop_shifted"};
			my $return=&parse_blast($out_tmp,$th_len,$start_e,$stop_e,$start_f,$seq_c,$start,$end,\%required);
			print "Return - $return\n";
			if ($return ne ""){
				my @t=sort {$$return{$b}{"id"} <=> $$return{$a}{"id"} || $$return{$b}{"length"} <=> $$return{$a}{"length"}} keys %{$return};
				print "LOOKING AT tRNA RETURNS\n";
				foreach my $hit (@t){
					print $hit."\n";
					if (!defined($store_att{$frag}{"tRNA"}) || ($store_att{$frag}{"tRNA"}{"id"}<$$return{$hit}{"id"}) || ($store_att{$frag}{"tRNA"}{"id"}==$$return{$hit}{"id"} && $$return{$hit}{"dr_length"} >$store_att{$frag}{"tRNA"}{"length"})){
						print "Att site in a tRNA - $hit - $$return{$hit}{pot_start}\n";
						$store_att{$frag}{"tRNA"}{"start"}=$$return{$hit}{"pot_start"};
						$store_att{$frag}{"tRNA"}{"end"}=$$return{$hit}{"pot_end"};
						$store_att{$frag}{"tRNA"}{"length"}=$$return{$hit}{"dr_length"};
						$store_att{$frag}{"tRNA"}{"id"}=$$return{$hit}{"id"};
					}
				}
			}
		}
		# Check CDS
		print "Now looking at integrase genes\n";
		# Look if we have integrase
		%putative_att=();
		for (my $i=$bottom_index-5;$i<=$top_index+5;$i++){
			if ($i<0 || $i>$#tab_genes){next;}
# 			print "$i - ".$store_gene{$frag}{$tab_genes[$i]}{"pfam"}."\n";
			if ($store_gene{$frag}{$tab_genes[$i]}{"pfam"}=~/integrase/i || $store_gene{$frag}{$tab_genes[$i]}{"pfam"} eq "rve" || $store_gene{$frag}{$tab_genes[$i]}{"pfam"} eq "Transposase_20"){
				$putative_att{$i}{"order"}=abs($mean-$i);
# 				print "$i is an integrase - $tab_genes[$i]\n";
			}
		}
		foreach my $index (sort {$putative_att{$b}{"order"} <=> $putative_att{$a}{"order"}} keys %putative_att){
			my $gene=$tab_genes[$index];
			print "Checking gene $gene\n";
			my $sens=$mean-$index;
			my $seq="";
			my $start_e;
			my $stop_e;
			if ($sens>0){# integrase is before center, we take the 1000 bp before the start
				$seq=substr($seq_c,$store_gene{$frag}{$gene}{"start_shifted"}-1000,1000);
				$start_e=$store_gene{$frag}{$gene}{"start_shifted"}-1000;
				if ($start_e<0){$start_e=0;}
				$stop_e=$store_gene{$frag}{$gene}{"start_shifted"};
			}
			else{ # integrase is after center, we take the 1000bp after the stop
				$seq=substr($seq_c,$store_gene{$frag}{$gene}{"stop_shifted"},1000);
				$start_e=$store_gene{$frag}{$gene}{"stop_shifted"};
				$stop_e=$store_gene{$frag}{$gene}{"stop_shifted"}+1000;
				if ($stop_e>length($seq_c)){$stop_e=$stop_e;}
			}
			open my $s1,">",$query_tmp;
			print $s1 ">Query\n$seq\n";
			close $s1;
			&run_cmd("blastn -query $query_tmp -db $db_tmp -out $out_tmp -task blastn-short -num_threads 1 -outfmt 6 -evalue 0.1");
			# We want hits with at least 90% ANI and an att site of 20bp
			my $return=&parse_blast($out_tmp,$th_len,$start_e,$stop_e,$start_f,$seq_c,$start,$end,\%required);
			if ($return ne ""){
				my @t=sort {$$return{$b}{"id"} <=> $$return{$a}{"id"} || $$return{$b}{"length"} <=> $$return{$a}{"length"}} keys %{$return};
				foreach my $hit (@t){
					if (!defined($store_att{$frag}{"Integrase"}) || ($store_att{$frag}{"Integrase"}{"id"}<$$return{$hit}{"id"}) || ($store_att{$frag}{"Integrase"}{"id"}==$$return{$hit}{"id"} && $$return{$hit}{"dr_length"} >$store_att{$frag}{"Integrase"}{"length"})){
						print "Att site with an integrase - $hit - $$return{$hit}{pot_start}\n";
						print "ID Is ".$$return{$hit}{"id"}."\n";
						$store_att{$frag}{"Integrase"}{"start"}=$$return{$hit}{"pot_start"};
						$store_att{$frag}{"Integrase"}{"end"}=$$return{$hit}{"pot_end"};
						$store_att{$frag}{"Integrase"}{"length"}=$$return{$hit}{"dr_length"};
						$store_att{$frag}{"Integrase"}{"id"}=$$return{$hit}{"id"};
					}
				}
			}
		}
		if ($store_att{$frag}{"tRNA"}{"length"}>=10 || $store_att{$frag}{"Integrase"}{"length"}>=10){
			$tag_canonical=1;
			my $new_start="";
			my $new_end="";
			my $type="";
			if ($store_att{$frag}{"tRNA"}{"length"}>=$store_att{$frag}{"Integrase"}{"length"}){
				print "## We found a canonical att site with a tRNA\n";
				print "New start / end in fragment are ".$store_att{$frag}{"tRNA"}{"start"}." -- ".$store_att{$frag}{"tRNA"}{"end"}."\n";
				$type="tRNA";
				$new_start=$store_att{$frag}{"tRNA"}{"start"}+$shift{$frag};
				$new_end=$store_att{$frag}{"tRNA"}{"end"}+$shift{$frag}-$store_att{$frag}{"tRNA"}{"length"}; # We don't take the repeat twice, i.e. we cut just before the att site
			}
			else{
				print "## We found a canonical att site with an integrase\n";
				print "New start / end in fragment are ".$store_att{$frag}{"Integrase"}{"start"}." -- ".$store_att{$frag}{"Integrase"}{"end"}."\n";
				$type="Integrase";
				$new_start=$store_att{$frag}{"Integrase"}{"start"}+$shift{$frag};
				$new_end=$store_att{$frag}{"Integrase"}{"end"}+$shift{$frag}-$store_att{$frag}{"Integrase"}{"length"}; # We don't take the repeat twice, i.e. we cut just before the att site
			}
			my $frag_length=$new_end-$new_start;
			print $frag." ".$start."-".$end." ".$frag_length." ".$check_frag{$frag}{"contig_length"}." -- ".$type."\n";
			print $s_out ">".$type.",".$frag."_".$new_start."-".$new_end.",".$frag_length.",".$check_frag{$frag}{"contig_length"}.",".$info_frag{$frag}{"confidence_score"}.",".$store_att{$frag}{$type}{"id"}.",".$store_att{$frag}{$type}{"length"}."\n";
			for (my $i=0;$i<=$#tab_genes;$i++){
				my $gene=$tab_genes[$i];
				if ($store_gene{$frag}{$gene}{"start"}>=$new_start && $store_gene{$frag}{$gene}{"start"}<=$new_end && $store_gene{$frag}{$gene}{"stop"}>=$new_start && $store_gene{$frag}{$gene}{"stop"}<=$new_end){
					print $s_out $store_gene{$frag}{$gene}{"line"}."\n";
				}
			}
		}
	}
	if ($tag_canonical==0){
		## No canonical att site found, look for "simple" Direct Repeats
		my $start_f=$start-5000;
		if ($start_f<0){$start_f=0;}
		print "We will take the fragment from $start_f to $start as database\n";
		my $seq_frag=substr($seq_c,$start_f,5000);
		open my $s1,">",$db_tmp;
		print $s1 ">Fragment\n$seq_frag\n";
		close $s1;
		&run_cmd("makeblastdb -in $db_tmp -title $db_tmp -dbtype nucl");
		my $stop_f=$end+5000;
		if ($stop_f>length($seq_c)){$stop_f=length($seq_c);}
		print "We will take the fragment from $end to $stop_f as query\n";
		my $seq_frag=substr($seq_c,$end,$stop_f-$end+1);
		open my $s1,">",$query_tmp;
		print $s1 ">Query\n$seq_frag\n";
		close $s1;
		print "Now looking at all repeats around the fragment\n";
		&run_cmd("blastn -query $query_tmp -db $db_tmp -out $out_tmp -task blastn-short -num_threads 1 -outfmt 6 -evalue 0.1");
		# We want hits with at least 95% ANI and an att site of 20bp
		my %hash=&parse_blast_untargeted($out_tmp,$th_len,$start_f,$end,$start);
		my @tab_pot=sort {$hash{$b}{"dr_length"} <=> $hash{$a}{"dr_length"} || $hash{$b}{"length"} <=> $hash{$a}{"length"}} keys %hash;
		if($#tab_pot<0){
			print "## no putative att sites identified\n";
			print $s_out $info_frag{$frag}{"data"}."\n";
		}
		else{
			my $tag=0;
			my $new_start=0;
			my $new_end="";
			my $hit_length="";
			my $hit_id="";
			my $type="";
			foreach my $hit (@tab_pot){
				if ($tag==2){last;}
				$tag=0;
# 				$type="DR";
				print "Potential att site $hit -> $hash{$hit}{pot_start}-$info_frag{$frag}{start}-$start====$end-$info_frag{$frag}{end}-$hash{$hit}{pot_end} --- $hash{$hit}{dr_length} / $hash{$hit}{id}\n";
				# Checking if this would be intergenic, intragenic, or add additional CDS - we will stop at the first one that is intergenic
				foreach my $gene (sort {$store_gene{$frag}{$a}{"start_shifted"} <=> $store_gene{$frag}{$b}{"start_shifted"}} keys %{$store_gene{$frag}}){
					if ($tag!=0){last;}
# 					print "$gene - $n --- ".$store_gene{$frag}{$gene}{"start_shifted"}." >= ".$hash{$hit}{"pot_start"}." && ".$store_gene{$frag}{$gene}{"stop_shifted"}." <= ".$hash{$hit}{"pot_end"}."\n";
					if ($store_gene{$frag}{$gene}{"start_shifted"}>=$start && $store_gene{$frag}{$gene}{"stop_shifted"}<=$end){
						print "\tGene $gene ".$store_gene{$frag}{$gene}{"start_shifted"}."-".$store_gene{$frag}{$gene}{"stop_shifted"}." already in the fragment\n";
					}
					elsif($store_gene{$frag}{$gene}{"start_shifted"}>=$hash{$hit}{"pot_start"} && $store_gene{$frag}{$gene}{"stop_shifted"}<=$hash{$hit}{"pot_end"}){
						print "\t##### Gene $gene ".$store_gene{$frag}{$gene}{"start_shifted"}."-".$store_gene{$frag}{$gene}{"stop_shifted"}." is an additional gene we just recruited ==== ".$store_gene{$frag}{$gene}{"type"}." --- ".$store_gene{$frag}{$gene}{"pfam"}." ".$store_gene{$frag}{$gene}{"pfam_score"}."\n";
						if ($store_gene{$frag}{$gene}{"type"} ne "CDS"){
							print "#### ==> hit $hit add an RNA, we remove\n";
							$tag=3;
						}
						elsif($store_gene{$frag}{$gene}{"pfam"} ne ""){
							if ($expected{$store_gene{$frag}{$gene}{"pfam"}}==1){}
							else{
								$tag=1;
								print "#### ==> hit $hit would add an unexpected PFAM, we remove\n";
							}
						}
					}
					elsif(($store_gene{$frag}{$gene}{"start_shifted"}>=$hash{$hit}{"pot_start"} && $store_gene{$frag}{$gene}{"start_shifted"}<=$hash{$hit}{"pot_end"}) || ($store_gene{$frag}{$gene}{"stop_shifted"}>=$hash{$hit}{"pot_start"} && $store_gene{$frag}{$gene}{"stop_shifted"}<=$hash{$hit}{"pot_end"})){
						print "\tGene $gene ".$store_gene{$frag}{$gene}{"start_shifted"}."-".$store_gene{$frag}{$gene}{"stop_shifted"}." overlap the potential ends\n";
						print "#### ==> hit $hit is in a CDS\n";
					}
				}
				if ($tag==0){
					$type="DR";
					$new_start=$hash{$hit}{"pot_start"}+$shift{$frag};
					$new_end=$hash{$hit}{"pot_end"}+$shift{$frag};
					$hit_length=$hash{$hit}{"dr_length"};
					$hit_id=$hash{$hit}{"id"};
					$tag=2;
					print "WE KEEP THIS ONE\n";
					# We found the one we want, no need to go further
				}
			}
			if ($new_start>0){
				print "##### New type $type\n";
				my $frag_length=$new_end-$new_start;
				print $frag." ".$start."-".$end." ".$frag_length." ".$check_frag{$frag}{"contig_length"}." -- ".$type."\n";
				print $s_out ">".$type.",".$frag."_".$new_start."-".$new_end.",".$frag_length.",".$check_frag{$frag}{"contig_length"}.",".$info_frag{$frag}{"confidence_score"}.",".$hit_id.",".$hit_length."\n";
				for (my $i=0;$i<=$#tab_genes;$i++){
					my $gene=$tab_genes[$i];
					if ($store_gene{$frag}{$gene}{"start"}>=$new_start && $store_gene{$frag}{$gene}{"start"}<=$new_end && $store_gene{$frag}{$gene}{"stop"}>=$new_start && $store_gene{$frag}{$gene}{"stop"}<=$new_end){
						print $s_out $store_gene{$frag}{$gene}{"line"}."\n";
					}
				}
			}
			else{
				print $s_out $info_frag{$frag}{"data"};
			}
		}
	}
}
close $s_out;

sub parse_blast{
	my $out_tmp=$_[0];
	my $th_len=$_[1];
	my $start_e=$_[2];
	my $stop_e=$_[3];
	my $decal=$_[4];
	my $seq_c=$_[5];
	my $start=$_[6];
	my $end=$_[7];
	my %required=%{$_[8]};
	my %return;
	open my $tsv,"<",$out_tmp;
	my $tag=0;
	my $i=0;
	while(<$tsv>){
		$i++;
		chomp($_);
		if ($tag==1){next;}
		my @tab=split("\t",$_);
		print "$_\n";
		my $match_len=$tab[7]-$tab[6]+1;
		if ($tab[2]>90 && $match_len>=$th_len){
			$tab[8]+=$decal;
			$tab[9]+=$decal;
			# if near the gene itself, it's a self hit
			if ($tab[9]<$tab[8]){my $t=$tab[9];$tab[9]=$tab[8];$tab[8]=$t;}
			if ($tab[8]==$start_e || $tab[9]==$stop_e || ($tab[8]>=$start_e && $tab[9]<=$stop_e)){
				print "$tab[8]-$tab[9] vs ".$start_e."-".$stop_e." == this is the gene/region vs itself\n";
			}
			else{
				# This is interesting => Look if this includes the fragment
				print "This $tab[8]-$tab[9] match is interesting (".$start_e."-".$stop_e.") -- id $tab[2]\n";
				my $id=$tab[2];
				my $pot_start;
				my $pot_extend;
				my $pot_stop;
				if ($tab[8]>$start_e){
					$pot_start=$start_e+$tab[6];
					$pot_extend=$start_e+$tab[7];
					$pot_stop=$tab[9];
				}
				else{
					$pot_start=$tab[8];
					$pot_extend=$tab[9];
					$pot_stop=$start_e+$tab[7];
				}
				my @tab=split("",$seq_c);
				my $extend=0;
				print "$pot_extend = $tab[$pot_extend] - $pot_stop = $tab[$pot_stop]\n";
				while($tab[$pot_extend] eq $tab[$pot_stop]){
# 					print "$pot_extend = $tab[$pot_extend] - $pot_stop = $tab[$pot_stop]\n";
					$extend++;
					$pot_extend++;
					$pot_stop++;
				}
				print "We extended it by $extend bp\n";
				my $diff_5=$start-$pot_extend;
				print "$pot_start => Compared to the predicted fragment ($start), there are $diff_5 additional/removed bp in 5'\n";
				my $diff_3=$pot_stop-$end;
				print "$pot_stop => Compared to the predicted fragment ($end), there are $diff_3 additional/removed bp in 3'\n";
				if ($diff_5>0 & $diff_3>0){
					print "We add bp, this is likely the correct one\n";
					$return{$i}{"pot_start"}=$pot_start;
					$return{$i}{"pot_end"}=$pot_stop;
					$return{$i}{"length"}=$pot_stop-$pot_start;
					$return{$i}{"id"}=$id;
					$return{$i}{"dr_length"}=$match_len+$extend;
					$return{$i}{"extend"}=$extend;
					$return{$i}{"diff_5"}=$diff_5;
					$return{$i}{"diff_3"}=$diff_3;
				}
				else{
					my $total_removed=0;
					if ($diff_5<0){$total_removed+=abs($diff_5);}
					if ($diff_3>0){$total_removed+=abs($diff_3);}
					if ($total_removed>0.5*($end-$start)){
						print "### We remove more than 50% of the fragment predicted, this is weird, likely not correct\n";
					}
					elsif($required{"start"}<$pot_start || $required{"stop"}<$pot_start || $required{"start"}>$pot_stop || $required{"stop"}>$pot_stop){
						print "### Pblm, this predicted prophage does not include the predicted Morph ($pot_start - $pot_stop vs $required{start} - $required{stop}\n";
					}
					else{
						print "This is possible \n";
						$return{$i}{"pot_start"}=$pot_start;
						$return{$i}{"pot_end"}=$pot_stop;
						$return{$i}{"length"}=$pot_stop-$pot_start;
						$return{$i}{"dr_length"}=$match_len+$extend;
						$return{$i}{"id"}=$id;
						$return{$i}{"extend"}=$extend;
						$return{$i}{"diff_5"}=$diff_5;
						$return{$i}{"diff_3"}=$diff_3;
					}
				}
			}
		}
	}
	close $tsv;
	return (\%return);
}

sub parse_blast_untargeted{
	my $out_tmp=$_[0];
	my $th_len=$_[1];
	my $start_f=$_[2];
	my $end=$_[3];
	my $start=$_[4];
	my %return;
	open my $tsv,"<",$out_tmp;
	my $tag=0;
	my $i=0;
	while(<$tsv>){
		chomp($_);
		if ($tag==1){next;}
		my @tab=split("\t",$_);
		print "### $_\n";
		my $match_len=$tab[7]-$tab[6]+1;
		if ($tab[2]>95 && $match_len>=$th_len){
			$tab[6]+=$end;
			$tab[7]+=$end;
			$tab[8]+=$start_f;
			$tab[9]+=$start_f;
			# This is interesting => Look if this includes the fragment
			print "This match between $tab[6]-$tab[7] and $tab[8]-$tab[9] is potentially interesting (".$start."-".$end.") Id -- $tab[2]\n";
			my $pot_start;
			my $pot_stop;
			$pot_start=$tab[8];
			if ($tab[9]<$pot_start){$pot_start=$tab[9];}
			$pot_stop=$tab[7];
			if ($tab[6] > $pot_stop){$pot_stop=$tab[6];}
			my $diff_5=$start-$pot_start;
			print "$pot_start => Compared to the predicted fragment, there are $diff_5 additional/removed bp in 5'\n";
			if ($diff_5<0){print "We remove some genes, so we don't take this one\n";next;}
			my $diff_3=$pot_stop-$end;
			if ($diff_3<0){print "We remove some genes, so we don't take this one\n";next;}
			print "$pot_stop => Compared to the predicted fragment, there are $diff_3 additional/removed bp in 3'\n";
			$return{$i}{"pot_start"}=$pot_start;
			$return{$i}{"pot_end"}=$pot_stop;
			$return{$i}{"dr_length"}=$tab[3];
			$return{$i}{"length"}=$pot_stop-$pot_start;
			$return{$i}{"id"}=$tab[2];
			$i++;
		}
	}
	close $tsv;
	return (%return);
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
