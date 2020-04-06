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
my $gb_file='';
my $path_pfam='';
my $n_cpu=2;
my $path_signalp='';
my $path_tmhmm='';
my $db_dir="Inovirus_db/";
GetOptions ('help' => \$h, 'h' => \$h, 'g=s'=>\$gb_file, 'p=s'=>\$path_pfam, 'd=s'=>\$db_dir, 't=s'=>\$n_cpu, , 'sp=s'=>\$path_signalp, 'th=s'=>\$path_tmhmm);
if ($h==1 || $gb_file eq "" || $path_pfam eq "" || $path_tmhmm eq "" || $path_signalp eq ""){ # If asked for help or did not set up any argument
	print "# Script to predict putative inovirus sequences from a gb file
#### Arguments : 
# -g : gb file of the genome to be analyzed
# -p : path to the pfam database hmm collection, i.e. Pfam-A.hmm (ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz)
# -sp : path to signal_p (http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp)
# -th : path to tmhmm (http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm)
#### Optional arguments
# -d : path to the Inovirus_db folder (default: Inovirus_db/)
# -t : number of threads used for hmmsearch and blast (default: 2)
#### Requirements:
# Bio::SeqIO
# Hmmsearch
# blastp
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
my $hmm_path=&run_cmd('which hmmsearch','quiet');
if ($hmm_path eq ""){die("pblm, we did not find the hmm path - $hmm_path, was the right module/conda environment loaded ?\n");}
my $blastp_path=&run_cmd('which blastp','quiet');
if ($blastp_path eq ""){die("pblm, we did not find the blastp path - $blastp_path, was the right module/conda environment loaded ?\n");}

## Take root name of the gb file
$gb_file=~/(.*)\.gb[^\.]*$/;
my $root_gb_file=$1;
if ($root_gb_file eq ""){die("Pblm with format file name $gb_file\n");}

print "### Compiling a fasta of all predicted proteins\n";
## make a fasta file of all proteins from the gb
my %prot_to_contig;
my $fasta_prots=$root_gb_file."_all_prots.faa";
if (-e $fasta_prots){
	print "Already a file $fasta_prots\n";
	## still have to load prot to contig
	my $io=Bio::SeqIO->new(-file => $gb_file, -format => 'genbank');
	while(my $seq=$io->next_seq){
		print "Taking proteins from ".$seq->id."\n";
		my @cds_features = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures;
		if (!defined($cds_features[0])){next;}
		my %seq_prot;
		if($cds_features[0]->has_tag('locus_tag')){
			%seq_prot = map { if ($_->has_tag('translation')){$_->get_tag_values('locus_tag') , $_->get_tag_values('translation')} else {$_->get_tag_values('locus_tag'), ""} } @cds_features;
		}
		elsif ($cds_features[0]->has_tag('protein_id')){
			%seq_prot = map { if ($_->has_tag('translation')){$_->get_tag_values('protein_id') , $_->get_tag_values('translation')} else {$_->get_tag_values('locus_tag'), ""} } @cds_features;
		}
		else{print "no protein id or locus tag, we are a little lost\n";}
		foreach my $id_prot (keys %seq_prot){
			$id_prot=~s/\s//g;
			$prot_to_contig{$id_prot}=$seq->id;
		}
		## Also linking tRNA and rRNA
		my @RNA_features = grep { $_->primary_tag eq 'rRNA' } $seq->get_SeqFeatures;
		push(@RNA_features,grep { $_->primary_tag eq 'tRNA' } $seq->get_SeqFeatures);
		if (!defined($RNA_features[0])){} # No rRNA, nothing to do
		else{
			foreach my $feat (@RNA_features){
				my $id_rna=($feat->get_tag_values('locus_tag'))[0];
				$prot_to_contig{$id_rna}=$seq->id;
# 				print "linking $id_rna to ".$seq->id."\n";
# 				<STDIN>;
			}
		}
	}
}
else{
	open my $s1,">",$fasta_prots;
	my $io=Bio::SeqIO->new(-file => $gb_file, -format => 'genbank');
	while(my $seq=$io->next_seq){
		print "Taking proteins from ".$seq->id."\n";
		my @cds_features = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures;
		if (!defined($cds_features[0])){next;}
		my %seq_prot;
		if($cds_features[0]->has_tag('locus_tag')){
			%seq_prot = map { if ($_->has_tag('translation')){$_->get_tag_values('locus_tag') , $_->get_tag_values('translation')} else {$_->get_tag_values('locus_tag'), ""} } @cds_features;
		}
		elsif ($cds_features[0]->has_tag('protein_id')){
			%seq_prot = map { if ($_->has_tag('translation')){$_->get_tag_values('protein_id') , $_->get_tag_values('translation')} else {$_->get_tag_values('locus_tag'), ""} } @cds_features;
		}
		else{
			print "no protein id or locus tag, we are a little lost\n";
		}
# 		
		foreach my $id_prot (keys %seq_prot){
			$id_prot=~s/\s//g;
			if ($seq_prot{$id_prot} eq ""){}
			else {
				print $s1 ">$id_prot\n$seq_prot{$id_prot}\n";
				print $id_prot."\n";
				$prot_to_contig{$id_prot}=$seq->id;
			}
		}
		my @RNA_features = grep { $_->primary_tag eq 'rRNA' } $seq->get_SeqFeatures;
		push(@RNA_features,grep { $_->primary_tag eq 'tRNA' } $seq->get_SeqFeatures);
		if (!defined($RNA_features[0])){} # No rRNA, nothing to do
		else{
			foreach my $feat (@RNA_features){
				my $id_rna=($feat->get_tag_values('locus_tag'))[0];
				$prot_to_contig{$id_rna}=$seq->id;
			}
		}
	}
	close $s1;
}

print "### Looking for pI-like proteins (inovirus marker gene) - hmmsearch against HMM profiles + blastp against singleton\n";
### Run Hmmsearch and blast 
my $out_hmm=$root_gb_file."_vs_Marker_db.tsv";
my $db_hmm=$db_dir."Final_marker_morph.hmm";
if (!(-e $out_hmm)){&run_cmd("hmmsearch -o /dev/null --noali --tblout $out_hmm --cpu $n_cpu $db_hmm $fasta_prots");}
else{print "$out_hmm already here\n";}
my $out_blast=$root_gb_file."_vs_ALV1.tsv";
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

if ($tag==0){die("no candidate morphogenesis found in this sequence\n");}

print "### Generating custom annotation file for each candidate fragment (i.e. 30kb around a pI-like protein) \n";
### Now going through the genomes, and extracting the fragments - we need for each candidate 30 kb around the candidate gene, we need the corresponding fasta files (nucleotide, and the set of proteins), we annotate these with PFAM, our custom database, and the coat predictor, and we generate a gff
my $io=Bio::SeqIO->new(-file => $gb_file, -format => 'genbank');
while(my $seq=$io->next_seq){
	my $seq_nucl=$seq->seq;
	print "Reading ".$seq->id."\n";
	my @cds_features = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures;
	if (!defined($cds_features[0])){next;}
	my %seq_prot;
	my %seq_start;
	my %seq_end;
	my %seq_strand;
	my %seq_product;
	my %seq_type;
	if($cds_features[0]->has_tag('locus_tag')){
		%seq_prot = map { if ($_->has_tag('translation')){$_->get_tag_values('locus_tag') , $_->get_tag_values('translation')} else {$_->get_tag_values('locus_tag'), ""} } @cds_features;
		%seq_start = map { $_->get_tag_values('locus_tag') , $_->start} @cds_features;
		%seq_end = map { $_->get_tag_values('locus_tag') , $_->end} @cds_features;
		%seq_product = map { if ($_->has_tag('product')){$_->get_tag_values('locus_tag') , $_->get_tag_values('product')} else {$_->get_tag_values('locus_tag'), ""} } @cds_features;
		%seq_strand = map { $_->get_tag_values('locus_tag') , $_->strand} @cds_features;
		%seq_type = map { $_->get_tag_values('locus_tag') , $_->primary_tag} @cds_features;
	}
	elsif ($cds_features[0]->has_tag('protein_id')){
		%seq_prot = map { if ($_->has_tag('translation')){$_->get_tag_values('protein_id') , $_->get_tag_values('translation')} else {$_->get_tag_values('protein_id'), ""} } @cds_features;
		%seq_start = map { $_->get_tag_values('protein_id') , $_->start} @cds_features;
		%seq_end = map { $_->get_tag_values('protein_id') , $_->end} @cds_features;
		%seq_product = map { if ($_->has_tag('product')){$_->get_tag_values('protein_id') , $_->get_tag_values('product')} else {$_->get_tag_values('protein_id'), ""} } @cds_features;
		%seq_strand = map { $_->get_tag_values('protein_id') , $_->strand} @cds_features;
		%seq_type = map { $_->get_tag_values('locus_tag') , $_->primary_tag} @cds_features;
	}
	else{
		print "no protein id or locus tag, we are a little lost\n";
	}
	## Also take rRNA and tRNA
	my @RNA_features = grep { $_->primary_tag eq 'rRNA' } $seq->get_SeqFeatures;
	push(@RNA_features,grep { $_->primary_tag eq 'tRNA' } $seq->get_SeqFeatures);
	if (!defined($RNA_features[0])){} # No rRNA, nothing to do
	else{
		print "We add some rRNA/tRNA\n";
# 		my %toto = map { $_->get_tag_values('locus_tag') , $_->start} @RNA_features;
# 		%seq_start = (%seq_start, %toto);
		%seq_start = (%seq_start, map { $_->get_tag_values('locus_tag') , $_->start} @RNA_features);
		%seq_end = (%seq_end, map { $_->get_tag_values('locus_tag') , $_->end} @RNA_features);
		%seq_strand = (%seq_strand, map { $_->get_tag_values('locus_tag') , $_->strand} @RNA_features);
		%seq_type = (%seq_type, map { $_->get_tag_values('locus_tag') , $_->primary_tag} @RNA_features);
		%seq_product = (%seq_product, map { if ($_->has_tag('product')){$_->get_tag_values('locus_tag') , $_->get_tag_values('product')} else {$_->get_tag_values('locus_tag'), ""} } @RNA_features);
	}
	
	foreach my $cand (keys %{$check_prot{$seq->id}}){
		my $b_s=$seq_start{$cand}-30000;
		if ($b_s<0){$b_s=0;}
		my $real_start=-1;
		my $b_e=$seq_end{$cand}+30000;
		my $real_end=-1;
		my $id_frag="frag_".$cand;
		my $out_file_fna=$root_gb_file."_".$id_frag."_nucl.fna";
		my $out_file_faa=$root_gb_file."_".$id_frag."_prots.faa";
		my $out_file_gff=$root_gb_file."_".$id_frag."_annot.gff";
		print "\t working with candidate $cand, $seq_start{$cand} / $seq_end{$cand} / $b_s / $b_e -> $out_file_fna / $out_file_faa / $out_file_gff\n";
		my @tab_selected;
		open my $s2,">",$out_file_faa;
		foreach my $prot (sort {$seq_start{$a} <=> $seq_start{$b}} keys %seq_start){
# 			print "Checking $prot - $seq_type{$prot} - $prot_to_contig{$prot} vs $prot_to_contig{$cand}\n";
			if ($prot_to_contig{$prot} eq $prot_to_contig{$cand}){
# 				print "Checking $prot - $seq_start{$prot}-$seq_end{$prot} vs $b_s-$b_e\n";
				if (($seq_start{$prot}>=$b_s && $seq_start{$prot}<=$b_e) || ($seq_end{$prot}>=$b_s && $seq_end{$prot}<=$b_e)){
					push(@tab_selected,$prot);
					print "\t\t We include feature $prot - $seq_type{$prot} - $seq_start{$prot} - $seq_end{$prot} - $real_start / $real_end\n";
					print $s2 ">".$prot."\n".$seq_prot{$prot}."\n";
					if ($real_start==-1){$real_start=$seq_start{$prot};}
					if ($seq_end{$prot}>$real_end){$real_end=$seq_end{$prot};}
				}
			}
		}
		close $s2;
		print "Real start = $real_start --- Real end = $real_end\n";
		open my $s1,">",$out_file_fna;
		print $s1 ">".$cand."\n";
		print $s1 substr($seq_nucl,$real_start,($real_end-$real_start+1));
		close $s1;
		## Do the hmm affi
		my $out_hmm_pfam=$root_gb_file."_".$id_frag."_vs_Pfam_db.tsv";
		if (!(-e $out_hmm_pfam)){&run_cmd("hmmsearch -o /dev/null --noali --tblout $out_hmm_pfam --cpu $n_cpu $path_pfam $out_file_faa");}
		else{print "$out_hmm_pfam already here\n";}
		my $out_hmm=$root_gb_file."_".$id_frag."_vs_InoSuperPCs.tsv";
		$db_hmm=$db_dir."Superclusters_db.hmm";
		if (!(-e $out_hmm)){&run_cmd("hmmsearch -o /dev/null --noali --tblout $out_hmm --cpu $n_cpu $db_hmm $out_file_faa");}
		else{print "$out_hmm already here\n";}
		my $out_blast=$root_gb_file."_".$id_frag."_vs_singletons.tsv";
		my $db_blast=$db_dir."Inoviridae_singletons";
		if (!(-e $out_blast)){&run_cmd("blastp -query $out_file_faa -db $db_blast -out $out_blast -num_threads $n_cpu -outfmt 6 -evalue 0.001");}
		else{print "$out_blast already here\n";}
		## Do the prediction of coat proteins
		my $out_coat_pred=$out_file_faa."_inovirus_coat_prediction.csv";
		if (!(-e $out_coat_pred)){&run_cmd($dirname."/Predict_inovirus_coat_proteins.pl -f $out_file_faa -sp $path_signalp -th $path_tmhmm","quiet");}
		else{print "$out_coat_pred already here\n";}
		if (!(-e $out_coat_pred)){die("Seems like there was a problem with Predict_inovirus_coat_proteins.pl -> we didn't get any output file \n");}
		## Load all annotation
		my %store_annot;
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
					$store_annot{$tab[0]}{"pfam"}{"name"}=join(" ",@tab[18..$#tab]);
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
		## Finalize the gff
		open my $s3,">",$out_file_gff;
		print $s3 "## ".$cand."\t".$real_start."\n";
		foreach my $prot (@tab_selected){
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
			if ($seq_strand{$prot}=~/^-/){$seq_strand{$prot}="-";}
			else{$seq_strand{$prot}="+";}
			if ($seq_type{$prot} eq "CDS"){
				print $s3 $seq->id."\tselected_feature\t".$seq_type{$prot}."\t".$seq_start{$prot}."\t".$seq_end{$prot}."\t.\t".$seq_strand{$prot}."\t0\tID=".$prot.";locus_tag=".$prot.";product=".$seq_product{$prot}.";pfam=".$code_pfam.";inoPC=".$code_ino.";coat_pred=".$code_coat."\n";
			}
			else{
				print $s3 $seq->id."\tselected_feature\t".$seq_type{$prot}."\t".$seq_start{$prot}."\t".$seq_end{$prot}."\t.\t".$seq_strand{$prot}."\t0\tID=".$prot.";locus_tag=".$prot.";product=".$seq_product{$prot}.";pfam=NA;inoPC=NA;coat_pred=NA\n";
			}
		}
		close $s3;
	}
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
