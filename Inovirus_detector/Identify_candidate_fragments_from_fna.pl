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
GetOptions ('help' => \$h, 'h' => \$h, 'f=s'=>\$original_fna_file, 'p=s'=>\$path_pfam, 'd=s'=>\$db_dir, 't=s'=>\$n_cpu, 'sp5=s'=>\$path_signalp5 , 'sp=s'=>\$path_signalp, 'th=s'=>\$path_tmhmm);
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
if (!($db_dir=~/\/$/)){$db_dir.="/";}
if (!(-d $db_dir)){die("pblm, we did not find the directory Inovirus_db -> $db_dir ?\n");}
my $dirname = dirname(__FILE__); ## Get the dir from which the script is called, so that we can get the full path of Predict_inovirus_coat_proteins.pl
# print "dir name: $dirname\n";<STDIN>;

## Test every program is here
my $prodigal_path=&run_cmd('which prodigal','quiet');
if ($prodigal_path eq ""){die("pblm, we did not find the prodigal path - $prodigal_path, was the right module/conda environment loaded ?\n");}
my $aragorn_path=&run_cmd('which aragorn','quiet');
if ($aragorn_path eq ""){die("pblm, we did not find the aragorn path - $aragorn_path, was the right module/conda environment loaded ?\n");}
my $hmm_path=&run_cmd('which hmmsearch','quiet');
if ($hmm_path eq ""){die("pblm, we did not find the hmm path - $hmm_path, was the right module/conda environment loaded ?\n");}
my $blastp_path=&run_cmd('which blastp','quiet');
if ($blastp_path eq ""){die("pblm, we did not find the blastp path - $blastp_path, was the right module/conda environment loaded ?\n");}

## Take root name of the input file
$original_fna_file=~/(.*)\.f[^\.]*$/;
my $root_fna_file=$1;
if ($root_fna_file eq ""){die("Pblm with format file name $original_fna_file\n");}
my $wdir="./";
if ($root_fna_file=~/(.*\/)[^\/]+/){$wdir=$1;}
print "Wdir: ".$wdir."\n";

## Clean sequence names
my $fna_file=$root_fna_file."_clean.fna";
if (-e $fna_file){print "Using existing $fna_file\n";}
else{
	open my $s1,">",$fna_file;
	open my $fna,"<",$original_fna_file;
	while(<$fna>){
		chomp($_);
		if ($_=~/^>(\S+)/){
			my $id=$1;
			$id=~s/[;:\/\.,\|\s\?!\*\%\+\@\=\{\}\~\[\]<>]/_/g;
			print $s1 ">".$id."\n";
		}
		else{
			print $s1 $_."\n";
		}
	}
	close $fna;
	close $s1;
}

print "### Compiling a fasta of all predicted proteins\n";
## make a fasta file of all proteins
my $fasta_prots=$root_fna_file."_all_cds.faa";
my $gff_prots=$root_fna_file."_all_cds.gff";
if (-e $fasta_prots){	print "Usint existing file $fasta_prots\n";}
else{
	### RUN PRODIGAL - GET FAA AND GFF
	### RUN INFERNAL - GET GFF (?)
	print "## Running prodigal -> $fasta_prots\n";
	&run_cmd("$prodigal_path -f gff -p meta -i $fna_file -a $fasta_prots -o $gff_prots");
}
### LOAD PROT TO CONTIG CLEANLY
my %prot_to_contig;
open my $faa,"<",$fasta_prots;
while(<$faa>){
	chomp($_);
	if ($_=~/^>(\S+)/){
		my $id=$1;
		if ($id=~/(.*)_\d+$/){
			$prot_to_contig{$id}=$1;
		}
            else{
                  print "## PBLM WITH $id format\n";
                  die("\n");
            }
	}
}
close $faa;


print "### Looking for pI-like proteins (inovirus marker gene) - hmmsearch against HMM profiles + blastp against singleton\n";
### Run Hmmsearch and blast
my $out_hmm=$root_fna_file."_vs_Marker_db.tsv";
my $db_hmm=$db_dir."Final_marker_morph.hmm";
if (!(-e $out_hmm)){&run_cmd("hmmsearch -o /dev/null --noali --tblout $out_hmm --cpu $n_cpu $db_hmm $fasta_prots");}
else{print "$out_hmm already here\n";}
my $out_blast=$root_fna_file."_vs_ALV1.tsv";
my $db_blast=$db_dir."Marker_ALV1";
if (!(-e $out_blast)){&run_cmd("blastp -query $fasta_prots -db $db_blast -out $out_blast -num_threads $n_cpu -outfmt 6 -evalue 0.001");}
else{print "$out_blast already here\n";}

### Look for candidates
my %check_prot;
my %store;
my $tag=0;
print "Reading $out_hmm ...\n";
open my $tsv,"<",$out_hmm;
while(<$tsv>){
	chomp($_);
	if ($_=~/^#/){next;}
	my @tab=split(" ",$_);
	if ($tab[5]>=30 && $tab[4]<=0.001){
		$check_prot{$prot_to_contig{$tab[0]}}{$tab[0]}=1;
		$tag=1;
		if (!defined($store{$tab[0]}) || $store{$tab[0]}{"score"}<$tab[5]){
			$store{$tab[0]}{"match"}=$tab[2];
			$store{$tab[0]}{"evalue"}=$tab[4];
			$store{$tab[0]}{"score"}=$tab[5];
			$store{$tab[0]}{"name"}=join(" ",@tab[18..$#tab]);
		}
	}
}
close $tsv;
print "Reading $out_blast ... \n";
open my $tsv,"<",$out_blast;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[11]>=50 && !defined($check_prot{$tab[0]})){
		$tag=1;
		$check_prot{$prot_to_contig{$tab[0]}}{$tab[0]}=1;
		$store{$tab[0]}{"match"}=$tab[1];
		$store{$tab[0]}{"evalue"}=$tab[10];
		$store{$tab[0]}{"score"}=$tab[11];
	}
}
close $tsv;

if ($tag==0){die("no candidate morphogenesis found in this input file\n");}

print "### Generating custom annotation file for each candidate fragment (i.e. 30kb around a pI-like protein) \n";
### Now going through the contigs, and extracting the fragments - we need for each candidate 30 kb around the candidate gene, we need the corresponding fasta files (nucleotide, and the set of proteins), we annotate these with PFAM, our custom database, and the coat predictor, and we generate a gff
my %store_annot;
my $c_contig="";
open my $gff,"<",$gff_prots;
while(<$gff>){
	chomp($_);
	if ($_=~/^#/){next;}
	my @tab=split("\t",$_);
	if ($tab[0] ne $c_contig && $c_contig ne ""){
		print "Processing $c_contig\n";
		foreach my $cand (keys %{$check_prot{$c_contig}}){
			### MAKE FUNCTIONS TO PROCESS
			&process($cand,$c_contig,\%store_annot);
		}
		### Now reset everything
		%store_annot=();
	}
	## NEED TO RECALCULATE ID
	my $id_prot=$tab[0];
	if ($tab[8]=~/ID=\d+_([^;]+)/){$id_prot.="_".$1;}
	else{die("Pblm with tab 8 $tab[8] in gff\n");}
	$store_annot{$id_prot}{"start"}=$tab[3];
	$store_annot{$id_prot}{"end"}=$tab[4];
	$store_annot{$id_prot}{"strand"}=$tab[6];
	$store_annot{$id_prot}{"type"}="CDS";
	$c_contig=$tab[0];
}
if ($c_contig ne ""){
	print "Processing $c_contig\n";
	foreach my $cand (keys %{$check_prot{$c_contig}}){
		### MAKE FUNCTIONS TO PROCESS
		print "\tGenerating file for candidate $cand (last contig)\n";
		&process($cand,$c_contig,\%store_annot);
	}
}


sub process{
	my $cand=$_[0];
	my $contig_id=$_[1];
	my %store_annot=%{$_[2]};
	if (!defined($store_annot{$cand}{"start"}) || (!defined($store_annot{$cand}{"end"}))){
		die("Pblm, did not find coordinates for candidate $cand\n");
	}
	my $b_s=$store_annot{$cand}{"start"}-30000;
	if ($b_s<0){$b_s=0;}
	my $real_start=-1;
	my $b_e=$store_annot{$cand}{"end"}+30000;
	my $real_end=-1;
	my $id_frag="frag_".$cand;
	my $out_file_fna=$root_fna_file."_".$id_frag."_nucl.fna";
	my $out_file_faa=$root_fna_file."_".$id_frag."_prots.faa";
	my $out_file_gff=$root_fna_file."_".$id_frag."_annot.gff";
	print "\t working with candidate $cand, ".$store_annot{$cand}{"start"}." / ".$store_annot{$cand}{"end"}." / $b_s / $b_e -> $out_file_fna / $out_file_faa / $out_file_gff\n";
	my @tab_selected;
	my %select_prot;
	foreach my $prot (sort {$store_annot{$a}{"start"} <=> $store_annot{$b}{"start"}} keys %store_annot){
		if ($prot_to_contig{$prot} eq $prot_to_contig{$cand}){
# 				print "Checking $prot - $seq_start{$prot}-$seq_end{$prot} vs $b_s-$b_e\n";
			if (($store_annot{$prot}{"start"}>=$b_s && $store_annot{$prot}{"start"}<=$b_e) || ($store_annot{$prot}{"end"}>=$b_s && $store_annot{$prot}{"end"}<=$b_e)){
				push(@tab_selected,$prot);
				$select_prot{$prot}=1;
				print "\t\t We include feature $prot - ".$store_annot{$prot}{"type"}." - ".$store_annot{$prot}{"start"}." - ".$store_annot{$prot}{"end"}." - $real_start / $real_end\n";
				if ($real_start==-1){$real_start=$store_annot{$prot}{"start"};}
				if ($store_annot{$prot}{"end"}>$real_end){$real_end=$store_annot{$prot}{"end"};}
			}
		}
	}
	# Fetching cds sequences
	my $tag=0;
	open my $faa,"<",$fasta_prots;
	open my $s2,">",$out_file_faa;
	while(<$faa>){
		chomp($_);
		if ($_=~/^>(\S+)/){
			my $id=$1;
			$tag=0;
			if ($select_prot{$id}==1){
				$tag=1;
				print $s2 ">".$id."\n";
			}
		}
		elsif($tag==1){
			print $s2 $_."\n";
		}
	}
	close $s2;
	close $faa;
	## Fragment fna
	print "Real start = $real_start --- Real end = $real_end\n";
	my $seq_nucl="";
	open my $fna,"<",$fna_file;
	$tag=0;
	while(<$fna>){
		chomp($_);
		if ($_=~/^>(\S+)/){
			my $id=$1;
			$tag=0;
			if ($id eq $contig_id){
				$tag=1;
			}

		}
		elsif($tag==1){
			$seq_nucl.=$_;
		}
	}
	close $fna;
	open my $s1,">",$out_file_fna;
	print $s1 ">".$cand."\n";
	print $s1 substr($seq_nucl,$real_start,($real_end-$real_start+1));
	close $s1;
	## Do the hmm affi
	my $out_hmm_pfam=$root_fna_file."_".$id_frag."_vs_Pfam_db.tsv";
	if (!(-e $out_hmm_pfam)){&run_cmd("hmmsearch -o /dev/null --noali --tblout $out_hmm_pfam --cpu $n_cpu $path_pfam $out_file_faa");}
	else{print "$out_hmm_pfam already here\n";}
	my $out_hmm=$root_fna_file."_".$id_frag."_vs_InoSuperPCs.tsv";
	$db_hmm=$db_dir."Superclusters_db.hmm";
	if (!(-e $out_hmm)){&run_cmd("hmmsearch -o /dev/null --noali --tblout $out_hmm --cpu $n_cpu $db_hmm $out_file_faa");}
	else{print "$out_hmm already here\n";}
	my $out_blast=$root_fna_file."_".$id_frag."_vs_singletons.tsv";
	my $db_blast=$db_dir."Inoviridae_singletons";
	if (!(-e $out_blast)){&run_cmd("blastp -query $out_file_faa -db $db_blast -out $out_blast -num_threads $n_cpu -outfmt 6 -evalue 0.001");}
	else{print "$out_blast already here\n";}
	## Do the prediction of coat proteins
	my $out_coat_pred=$out_file_faa."_inovirus_coat_prediction.csv";
	if (!(-e $out_coat_pred)){
		if ($path_signalp5 ne ""){&run_cmd($dirname."/Predict_inovirus_coat_proteins.pl -f $out_file_faa -sp5 $path_signalp5 -th $path_tmhmm -w $wdir","quiet");}
		else{&run_cmd($dirname."/Predict_inovirus_coat_proteins.pl -f $out_file_faa -sp $path_signalp -th $path_tmhmm -w $wdir","quiet");}
	}
	else{print "$out_coat_pred already here\n";}
	if (!(-e $out_coat_pred)){die("Seems like there was a problem with Predict_inovirus_coat_proteins.pl -> we didn't get any output file \n");}
	## Do the prediction of tRNA with Infernal
	my $out_aragorn=$root_fna_file."_".$id_frag."_tRNAs.tsv";
	if (!(-e $out_aragorn)){&run_cmd("$aragorn_path -t $out_file_fna -o $out_aragorn -w","quiet");}
	else{print "$out_aragorn already here\n";}
	## Load all annotation
	print "Reading $out_hmm_pfam ... \n";
	open my $tsv,"<",$out_hmm_pfam;
	while(<$tsv>){
		chomp($_);
		if ($_=~/^#/){next;}
		my @tab=split(" ",$_);
		if ($tab[5]>=30 && $tab[4]<=0.001){
			if (!defined($store_annot{$tab[0]}{"pfam"}) || $store_annot{$tab[0]}{"pfam"}{"score"}<$tab[5]){
				$store_annot{$tab[0]}{"pfam"}{"match"}=$tab[2];
				$store_annot{$tab[0]}{"pfam"}{"evalue"}=$tab[4];
				$store_annot{$tab[0]}{"pfam"}{"score"}=$tab[5];
				$store_annot{$tab[0]}{"product"}="Putative ".$tab[3];
			}
		}
	}
	close $tsv;

	print "Reading $out_hmm ... \n";
	open my $tsv,"<",$out_hmm;
	while(<$tsv>){
		chomp($_);
		if ($_=~/^#/){next;}
		my @tab=split(" ",$_);
		if ($tab[5]>=30 && $tab[4]<=0.001){
			if (!defined($store_annot{$tab[0]}{"inoPC"}) || $store_annot{$tab[0]}{"inoPC"}{"score"}<$tab[5]){
				$store_annot{$tab[0]}{"inoPC"}{"match"}=$tab[2];
				$store_annot{$tab[0]}{"inoPC"}{"evalue"}=$tab[4];
				$store_annot{$tab[0]}{"inoPC"}{"score"}=$tab[5];
				$store_annot{$tab[0]}{"inoPC"}{"name"}=join(" ",@tab[18..$#tab]);
			}
		}
	}
	close $tsv;


	print "Reading $out_blast ... \n";
	open my $tsv,"<",$out_blast;
	while(<$tsv>){
		chomp($_);
		my @tab=split("\t",$_);
		if ($tab[11]>=50){
			if (!defined($store_annot{$tab[0]}{"inoSglton"}) || $store_annot{$tab[0]}{"inoSglton"}{"score"}<$tab[11]){
				$store_annot{$tab[0]}{"inoSglton"}{"match"}=$tab[1];
				$store_annot{$tab[0]}{"inoSglton"}{"evalue"}=$tab[10];
				$store_annot{$tab[0]}{"inoSglton"}{"score"}=$tab[11];
			}
		}
	}
	close $tsv;

	print "Reading $out_coat_pred .. \n";
	open my $csv,"<",$out_coat_pred;
	while(<$csv>){
		chomp($_);
		if ($_=~/^#/){next;}
		my @tab=split(",",$_);
		$store_annot{$tab[0]}{"coat_pred"}=$tab[1];
	}
	close $csv;

	print "Reading $out_aragorn .. \n";
	open my $tsv,"<",$out_aragorn;
	while(<$tsv>){
		chomp($_);
		my @tab=split(" ",$_);
		if ($tab[1]=~/^tRNA/){
			my $id=$tab[1]."_".$tab[0]; ## tab[0] has the detection #, so it will be a uniq ID every time
			push(@tab_selected,$id);
			$store_annot{$id}{"type"}="tRNA";
			$store_annot{$id}{"product"}=$tab[1];
			if ($tab[2]=~/^\[(\d+),(\d+)\]/){
				$store_annot{$id}{"start"}=$real_start+$1; ## Has to shift all coordinates
				$store_annot{$id}{"end"}=$real_start+$2;
				$store_annot{$id}{"strand"}="+";
			}
			elsif($tab[2]=~/^c\[(\d+),(\d+)\]/){
				$store_annot{$id}{"start"}=$real_start+$1;
				$store_annot{$id}{"end"}=$real_start+$2;
				$store_annot{$id}{"strand"}="-";
			}
			else{
				print "PBLM ARAGORN $_\n";
				print $tab[2]."\n";
				die("\n");
			}
		}
	}
	close $tsv;


	## Finalize the gff
	open my $s3,">",$out_file_gff;
	print $s3 "## ".$cand."\t".$real_start."\n";
	foreach my $prot (sort {$store_annot{$a}{"start"} <=> $store_annot{$b}{"start"}} @tab_selected){
		my $code_pfam="NA";
		if (defined($store_annot{$prot}{"pfam"}{"match"})){
			$code_pfam=$store_annot{$prot}{"pfam"}{"match"}." ".$store_annot{$prot}{"pfam"}{"score"};
		}
		my $code_ino="NA";
		if (defined($store_annot{$prot}{"inoPC"}{"match"})){
			$code_ino=$store_annot{$prot}{"inoPC"}{"match"}." ".$store_annot{$prot}{"inoPC"}{"score"};
		}
		elsif (defined($store_annot{$prot}{"inoSglton"}{"match"})){
			$code_ino=$store_annot{$prot}{"inoSglton"}{"match"}." ".$store_annot{$prot}{"inoSglton"}{"score"};
		}
		my $code_coat=$store_annot{$prot}{"coat_pred"};
		if ($store_annot{$prot}{"strand"}=~/^-/){$store_annot{$prot}{"strand"}="-";}
		else{$store_annot{$prot}{"strand"}="+";}
		if ($store_annot{$prot}{"type"} eq "CDS"){
			print $s3 $contig_id."\tselected_feature\t".$store_annot{$prot}{"type"}."\t".$store_annot{$prot}{"start"}."\t".$store_annot{$prot}{"end"}."\t.\t".$store_annot{$prot}{"strand"}.
			"\t0\tID=".$prot.";locus_tag=".$prot.";product=".$store_annot{$prot}{"product"}.";pfam=".$code_pfam.";inoPC=".$code_ino.";coat_pred=".$code_coat."\n";
		}
		else{
			print $s3 $contig_id."\tselected_feature\t".$store_annot{$prot}{"type"}."\t".$store_annot{$prot}{"start"}."\t".$store_annot{$prot}{"end"}."\t.\t".$store_annot{$prot}{"strand"}.
			"\t0\tID=".$prot.";locus_tag=".$prot.";product=".$store_annot{$prot}{"product"}.";pfam=NA;inoPC=NA;coat_pred=NA\n";
		}
	}
	close $s3;
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
		chomp($out);
		return($out);
	}
	else{
		print " ### dummy run\n";
	}
}
