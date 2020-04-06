#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h='';
my $dir='';
my $wdir='';
my $n_cpu=2;
my $input_file="";
my $db_dir="";
my $out_dir="";
my $tag_tmhmm=1;
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$input_file, 'd=s'=>\$db_dir, 'o=s'=>\$out_dir);
if ($h==1 || $input_file eq "" || $db_dir eq "" || $out_dir eq ""){ # If asked for help or did not set up any argument
	print "# Script to classify a set of Inovirus genomes
# -i: fasta file of inovirus genomes to be classified
# -d: path to the \"Ino_classifier_db/\" database
# -o: output directory (will be created if needed)
";
	die "\n";
}
my $prodigal_path=&run_cmd('which prodigal','veryquiet');
my $hmm_path=&run_cmd('which hmmsearch','veryquiet');
my $blast_path=&run_cmd('which blastp','veryquiet');

if (!(-d $out_dir)){&run_cmd("mkdir $out_dir");}

## List genomes
my %info_genome;
open my $fa,"<",$input_file;
while(<$fa>){
      chomp($_);
      if ($_=~/^>(\S+)/){
            my $id=$1;
            $info_genome{$id}{"selected"}=1;
      }
}
close $fa;

# Step 1: prodigal protein prediction
my $id_run="default_id";
if ($input_file=~/([^\/]+)\.f[^\/\s]+$/){$id_run=$1;}
my $faa_file=$out_dir.$id_run."_cds.faa";
if (-e $faa_file){print "Use existing $faa_file\n";}
else{
      print "## Running prodigal -> $faa_file\n";
      &run_cmd("$prodigal_path -p meta -i $input_file -a $faa_file")
}

# Step 2: hmmsearch against custom Morph references
my $db_morph=$db_dir."/pI_PCs_db_annot.hmm";
my $hmm_morph=$out_dir.$id_run."_vs_morph_hmm.tab";
my $hmm_morph_by_dom=$out_dir.$id_run."_vs_morph_hmm_domain.tab";
if (-e $hmm_morph){print "Use existing $hmm_morph\n";}
else{
      print "## Running hmmsearch -> $hmm_morph ($hmm_path)\n";
      &run_cmd("$hmm_path --tblout $hmm_morph --domtblout $hmm_morph_by_dom --cpu $n_cpu -o /dev/null --noali $db_morph $faa_file");
}
my %store_morph;
my %check_morph;
open my $tsv,"<",$hmm_morph;
while(<$tsv>){
	chomp($_);
	if ($_=~/^#/){next;}
	my @tab=split(" ",$_);
      if ($tab[2] eq "Zot"){next;} ## Zot is used only to detect ATPase domain boundary, so actually not used here
	if ($tab[5]>=50 && $tab[8]>=50){
            if ($tab[0]=~/(.*)_\d+$/){
                  my $genome=$1;
      		if (!defined($store_morph{$genome}{$tab[0]}) || $store_morph{$genome}{$tab[0]}{"score"}<$tab[5]){
      			$store_morph{$genome}{$tab[0]}{"hit"}=$tab[2];
      			$store_morph{$genome}{$tab[0]}{"evalue"}=$tab[4];
      			$store_morph{$genome}{$tab[0]}{"score"}=$tab[5];
      		}
                  $check_morph{$tab[0]}=1;
            }
            else{
                  print "## PBLM WITH $tab[0] format\n";
                  die("\n");
            }
	}
}
close $tsv;

# Step 3: blastp against proteins from reference genomes
my %store_blast;
my %store_length;
my $db_prot=$db_dir."/Ref_database";
my $info_file=$db_dir."/Info_refs_wMorph.tsv";
my $blast_result=$out_dir.$id_run."_vs_Ino_refs.tab";
if (-e $blast_result){print "Use existing $blast_result\n";}
else{
      print "## Running blastp -> $blast_result ($blast_path)\n";
	my $test_db=$db_prot.".psq";
	if (!(-e $test_db)){&run_cmd("makeblastdb -in $db_prot.faa -out $db_prot -dbtype prot");}
      &run_cmd("$blast_path -query $faa_file -db $db_prot -out $blast_result -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident positive\" -evalue 0.001 -num_threads $n_cpu");
}
my %affi_ref;
my %check_morph_ref;
open my $tsv,"<",$info_file;
while(<$tsv>){
      chomp($_);
      if ($_=~/^#/){next;}
      my @tab=split("\t",$_);
      $affi_ref{$tab[0]}{"name"}=$tab[1];
      $affi_ref{$tab[0]}{"family"}=$tab[2];
      $affi_ref{$tab[0]}{"genus"}=$tab[3];
	if ($tab[4] ne "NA"){$check_morph_ref{$tab[4]}=1;}
}
close $tsv;

open my $tsv,"<",$blast_result;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[10]<=0.001){ ## Maximum e-value => 0.001 ### MAY WANT TO CHANGE TO A SCORE TH ONE DAY ?? BUT MAYBE NOT IF FIXED DATABASE
            if ($tab[0]=~/(.*)_\d+$/){
                  my $g1=$1;
                  my $p1=$tab[0];
                  my @t2=split(/\|/,$tab[1]);
                  my $p2=$t2[0];
                  my $g2=$t2[1];
			### WE KNOW WE WILL HAVE SPURIOUS HITS ON pI, SO WE ADD A THRESHOLD ON SCORE FOR THIS ONE
			if (($check_morph{$tab[0]}==1 || $check_morph_ref{$tab[1]}==1) && $tab[11]<50){
				# print "We don't consider hit between $tab[0] and $tab[1] because it includes Morph proteins and has a score < 50 ($tab[11])\n";
				next;
			}
                  if (defined($store_blast{$g1}{$g2}{"seen"}{$p1})){ ## We already had a hit for this protein in this genome
                  }
                  elsif($store_blast{$g1}{$g2}{"seen"}{$p2}){## We already had a hit for this protein in this genome
                  }
                  else{ ## This is a new pair
                        $store_blast{$g1}{$g2}{"seen"}{$p1}=1;
                        $store_blast{$g1}{$g2}{"seen"}{$p2}=1;
                        $store_blast{$g1}{$g2}{"total"}++;
                  }
      		$store_blast{$g1}{$g2}{"match"}{$p1}{$p2}{"ident"}+=$tab[14];
      		if (!defined($store_length{$p1})){
      			$store_length{$p1}=$tab[12];
      		}
      		if (!defined($store_length{$p2})){
      			$store_length{$p2}=$tab[13];
      		}
            }
	      else{
                  print "## PBLM WITH $tab[0] format\n";
                  die("\n");
            }
      }
}
close $tsv;

# Final: Do an affiliation based on best hit of Morph + 2-gene affi to any refernce genome
my $max=5;
my $final_out=$out_dir.$id_run."_final_classification.tsv";
my $final_out_bis=$out_dir.$id_run."_final_classification_forR.tsv";
open my $s_final,">",$final_out;
print $s_final "## Genome\tClosest reference\tCorresponding family\tCorresponding genus\tList of significant hits\n";
foreach my $genome (sort keys %info_genome){
      ## 2 genes affiliation
      my @list_hits = sort {$store_blast{$genome}{$b}{"total"} <=> $store_blast{$genome}{$a}{"total"}} keys %{$store_blast{$genome}};
      my %tmp;
      for (my $i=0;$i<$max;$i++){
            my $hit=$list_hits[$i];
            my $total=$store_blast{$genome}{$hit}{"total"};
            ## Calculate AAI
            my $cumulated_aai=0;
            foreach my $prot (sort keys %{$store_blast{$genome}{$hit}{"match"}}){
                  my @prot_hit = sort { $store_blast{$genome}{$hit}{"match"}{$prot}{$b}{"ident"} <=> $store_blast{$genome}{$hit}{"match"}{$prot}{$a}{"ident"} } keys %{$store_blast{$genome}{$hit}{"match"}{$prot}};
                  my $aai=$store_blast{$genome}{$hit}{"match"}{$prot}{$prot_hit[0]}{"ident"}/$store_length{$prot}*100;
                  $cumulated_aai+=$aai;
            }
            if ($total>=2){
                  $tmp{$hit}=$cumulated_aai;
            }
            if ($i==0 && $total<2){
                  $tmp{"new_family"}=1;
            }
      }
      my $hit_blast="NA";
      my $affi_blast_family="NA";
      my $affi_blast_genus="NA";
      my @t=sort {$tmp{$b} <=> $tmp{$a}} keys %tmp;
      if (!defined($t[0]) || $t[0] eq "new_family"){
            $affi_blast_family="Unclassified";
		$affi_blast_genus="Unclassified";
      }
	elsif($#t==0){ ## Only one reference with 2 hits, we take
		$hit_blast=$t[0];
		$affi_blast_family=$affi_ref{$t[0]}{"family"};
		$affi_blast_genus=$affi_ref{$t[0]}{"genus"};
	}
      elsif ($#t>0){
            $hit_blast=$t[0];
            $affi_blast_family=$affi_ref{$t[0]}{"family"};
            $affi_blast_genus=$affi_ref{$t[0]}{"genus"};
		for (my $i=1;$i<=$#t;$i++){
			if ($affi_ref{$t[$i]}{"family"} eq $affi_blast_family){} ## All good, same family
			else{
				if (($tmp{$t[$i]}/$tmp{$hit_blast})>0.9){
					## Near-identical hit, we can't conclude anything
			            $affi_blast_family="Unclassified";
					$affi_blast_genus="Unclassified";
					$i=$#t+1; ## We can end there
				}
			}
		}
      }
	my $hit_list="";
	for (my $i=0;$i<5;$i++){
		if (defined($t[$i])){
			$hit_list.=$t[$i]." ".$tmp{$t[$i]}.";";
		}
	}
	chop($hit_list);
      print $s_final $genome."\t".$hit_blast."\t".$affi_blast_family."\t".$affi_blast_genus."\t".$hit_list."\n";
      print $genome."\t".$hit_blast."\t".$affi_blast_family."\t".$affi_blast_genus."\t".$hit_list."\n";
}
close $s_final;
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
