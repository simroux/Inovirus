#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h='';
my $cmd='';
my $out='';
my $path_signalp='/global/u2/s/simroux/Utils/signalp/signalp-4.1/signalp';
my $path_signalp5='/global/u2/s/simroux/Utils/signalp/signalp-5.0/signalp-5.0b/';
my $path_tmhmm='/global/u2/s/simroux/Utils/tmhmm-2.0c/bin/tmhmm';
my $fa_file='';
my $wdir='./';
GetOptions ('help' => \$h, 'h' => \$h, 'f=s'=>\$fa_file, 'sp=s'=>\$path_signalp, 'sp5=s'=>\$path_signalp5, 'th=s'=>\$path_tmhmm, 'w=s'=>\$wdir);
if ($h==1 || $fa_file eq ""  || $path_tmhmm eq "" || ($path_signalp eq "" && $path_signalp5 eq "")){ # If asked for help or did not set up any argument
	print "# Script to predict putative inovirus coat proteins
#### Arguments :
# -f : fasta file of the proteins to be analyzed
# -th : path to tmhmm (http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm)
# -sp : path to signal_p <v4 (http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp)
# or
# -sp5 : path to signal_p v5 folder (www.cbs.dtu.dk/services/SignalP/portable.php)
#### Requirements:
# Tmhmm
# SignalP
####
";
	die "\n";
}


if (!($wdir=~/\/$/)){$wdir.="/";}

my $tmp_dir_tmhmm=$wdir."TmpdirTMHMM";
if (!(-d $tmp_dir_tmhmm)){`mkdir $tmp_dir_tmhmm`;}
my $out_file=$fa_file."_inovirus_coat_prediction.csv";

open my $fa,"<",$fa_file;
my @tab_genes;
while(<$fa>){
	chomp($_);
	if ($_=~/^>(\S+)/){push(@tab_genes,$1);}
}
close $fa;


print "Run Tmhmm on native proteins\n";
my $tmhmm=$wdir."out_tmp_tmhmm.txt";
&run_tmhmm($fa_file,$tmhmm);
my $c_c="";
my %info;
open my $txt,"<",$tmhmm;
while(<$txt>){
	chomp($_);
	if ($_=~/^# (\S+)\s+Length:\s+(\d+)/){
		$c_c=$1;
		my $length=$2;
# 		my @t=split(/\./,$c_c);
		$info{$c_c}{"length"}=$length
	}
	elsif($_=~/^# (\S+)\s+Number of predicted TMHs:\s+(\d+)/){
		$c_c=$1;
		my $tm=$2;
# 		my @t=split(/\./,$c_c);
		$info{$c_c}{"tm"}=$tm;
	}
}
close $txt;

my %store_signal;
my $tmp_matured=$wdir."Tmpmatured.faa";
my $result_gp=$wdir."out_tmp_signalp_gramp";
my $matured_gp=$wdir."Tmpmatured_gp.faa";
my $result_gn=$wdir."out_tmp_signalp_gramn";
my $matured_gn=$wdir."Tmpmatured_gn.faa";
if ($path_signalp5 eq ""){
	print "Generating matured proteins\n";
	print "Running SignalP with gram+ model\n";
	&run_signalp($fa_file,"gram+",$result_gp,$matured_gp);
	open my $tsv,"<",$result_gp;
	my $test=0;
	while(<$tsv>){
		chomp($_);
		if ($_ ne ""){$test++;}
		if ($_=~/^#/){next;}
		my @tab=split(" ",$_);
		$store_signal{$tab[0]}{"gp_score"}=$tab[3];
	}
	close $tsv;
	if ($test<1){
		die("Pblm with signalp - no result ? \n");
	}
	print "Running SignalP with gram- model\n";
	&run_signalp($fa_file,"gram-",$result_gn,$matured_gn);
	open my $tsv,"<",$result_gn;
	$test=0;
	while(<$tsv>){
		chomp($_);
		if ($_ ne ""){$test++;}
		if ($_=~/^#/){next;}
		my @tab=split(" ",$_);
		$store_signal{$tab[0]}{"gn_score"}=$tab[3];
	}
	close $tsv;
	if ($test<1){
		die("Pblm with signalp - no result for gram - ? \n");
	}
}
else{
	print "Running SignalP5 with gram+ model\n";
	&run_signalp5($fa_file,"gram+",$result_gp,$matured_gp);
	my %store_signal;
	open my $tsv,"<",$result_gp;
	my $test=0;
	while(<$tsv>){
		chomp($_);
		if ($_ ne ""){$test++;}
		if ($_=~/^#/){next;}
		my @tab=split(" ",$_);
		$store_signal{$tab[0]}{"gp_score"}=$tab[2];
	}
	close $tsv;
	if ($test<1){
		die("Pblm with signalp - no result ? \n");
	}
	print "Running SignalP5 with gram- model\n";
	if (-e $result_gn){}
	else{&run_signalp5($fa_file,"gram-",$result_gn,$matured_gn);}
	open my $tsv,"<",$result_gn;
	$test=0;
	while(<$tsv>){
		chomp($_);
		if ($_ ne ""){$test++;}
		if ($_=~/^#/){next;}
		my @tab=split(" ",$_);
		$store_signal{$tab[0]}{"gn_score"}=$tab[2];
	}
	close $tsv;
}

my %store_seq;
my $tag=0;
open my $s1,">",$tmp_matured;
open my $fa,"<",$matured_gp;
while(<$fa>){
	chomp($_);
	if ($_=~/^>(.*)/){
		$tag=0;
		my $id=$1;
		if ($path_signalp5 eq ''){
			if ($id=~/^(\S+).*; MatureChain: (\d+-\d+)/){
				my $c_c=$1;
				my $chain=$2;
				if (!defined($store_signal{$c_c}{"gn_score"}) || $store_signal{$c_c}{"gp_score"}>=$store_signal{$c_c}{"gn_score"}){
					$tag=1;
					print $s1 ">".$c_c." gramp,".$chain."\n";
				}
			}
			else{
				print "# ___ Pblm with $_\n";
				<STDIN>;
			}
		}
		else{
			if (!defined($store_signal{$id}{"gn_score"}) || $store_signal{$id}{"gp_score"}>=$store_signal{$id}{"gn_score"}){
				$tag=1;
				print $s1 ">".$c_c." gramp\n";
			}
		}
	}
	elsif ($tag==1){
		print $s1 "$_\n";
	}
}
close $fa;
open my $fa,"<",$matured_gn;
while(<$fa>){
	chomp($_);
	if ($_=~/^>(.*)/){
		$tag=0;
		my $id=$1;
		if ($path_signalp5 eq ''){
			if ($id=~/^(\S+).*; MatureChain: (\d+-\d+)/){
				my $c_c=$1;
				my $chain=$2;
				if (!defined($store_signal{$c_c}{"gp_score"}) || $store_signal{$c_c}{"gp_score"}<$store_signal{$c_c}{"gn_score"}){
					$tag=1;
					print $s1 ">".$c_c." gramn,".$chain."\n";
				}
			}
			else{
				print "# ___ Pblm with $_\n";
				<STDIN>;
			}
		}
		else{
			if (!defined($store_signal{$id}{"gp_score"}) || $store_signal{$id}{"gp_score"}<=$store_signal{$id}{"gn_score"}){
				$tag=1;
				print $s1 ">".$c_c." gramn\n";
			}
		}
	}
	elsif ($tag==1){
		print $s1 "$_\n";
	}
}
close $fa;
close $s1;

print "run TMHMM on matured proteins\n";
my $tmhmm_2=$wdir."out_tmp_tmhmm_matured.txt";
&run_tmhmm($tmp_matured,$tmhmm_2);
my %info_2;
open my $txt,"<",$tmhmm_2;
while(<$txt>){
	chomp($_);
	if ($_=~/^# (\S+)\s+Length:\s+(\d+)/){
		$c_c=$1;
		my $length=$2;
# 		my @t=split(/\./,$c_c);
		$info_2{$c_c}{"length"}=$length
	}
	elsif($_=~/^# (\S+)\s+Number of predicted TMHs:\s+(\d+)/){
		$c_c=$1;
		my $tm=$2;
# 		my @t=split(/\./,$c_c);
		$info_2{$c_c}{"tm"}=$tm;
	}
}
close $txt;
my %pred;
foreach my $prot (@tab_genes){
	if ($info_2{$prot}{"length"}>30 && $info{$prot}{"length"}<60 && $info_2{$prot}{"tm"}==1){
		print "$prot\t".$info_2{$prot}{"length"}."\t".$info_2{$prot}{"tm"}."\tPerfect after SignalP\n";
		$pred{$prot}="Perfect_sp";
	}
	elsif($info{$prot}{"length"}>30 && $info{$prot}{"length"}<60 && $info{$prot}{"tm"}==1){
		print "$prot\t".$info{$prot}{"length"}."\t".$info{$prot}{"tm"}."\tPerfect natively\n";
		$pred{$prot}="Perfect";
	}
	elsif ($info{$prot}{"length"}>30 && $info{$prot}{"length"}<120 && $info{$prot}{"tm"}>=1 && $info{$prot}{"tm"}<=2){
		print "$prot\t".$info{$prot}{"length"}."\t".$info{$prot}{"tm"}."\tSome signal but not perfect\n";
		$pred{$prot}="Signal";
	}
}
open my $final_out,">",$out_file;
print $final_out "### Perfect: length of 30 to 60 aa, with a single transmembrane domain, likely an inovirus coat protein\n";
print $final_out "### Perfect_sp: length of 30 to 60 aa, with a single transmembrane domain, likely an inovirus coat protein\n";
print $final_out "### Signal: length of 30 to 120 aa, with a single transmembrane domain or two transmembrane domains, possibly an inovirus coat protein\n";
print $final_out "### NA: most likely not an inovirus coat protein\n";
print $final_out "Protein,Inovirus coat prediction\n";
foreach my $gene (@tab_genes){
	if (!defined($pred{$gene})){$pred{$gene}="NA";}
	print $final_out $gene.",".$pred{$gene}."\n";
}
close $final_out;

# print "### Cleaning up \n";
# &run_cmd("rm $tmhmm $tmhmm_2 $tmp_matured $result_gp $matured_gp $result_gn $matured_gn");




sub run_tmhmm{
	my $in_file=$_[0];
	my $out_file=$_[1];
	my $tmp_dir=$tmp_dir_tmhmm;
	if (!(-d $tmp_dir)){&run_cmd("mkdir $tmp_dir");}
	&run_cmd("echo '' > $out_file");
	open my $fa,"<",$in_file;
	my $c_c="";
	while(<$fa>){
		chomp($_);
		if ($_=~/^>/){
			if ($c_c ne ""){
				&run_cmd("echo \"$c_c\" | $path_tmhmm -workdir $tmp_dir >> $out_file","quiet");
				&run_cmd("rm -rf $tmp_dir/*.*");
			}
			$c_c=$_."\n";
		}
		else{
			$c_c.=$_."\n";
		}
	}
	close $fa;
	&run_cmd("echo \"$c_c\" | $path_tmhmm -workdir $tmp_dir >> $out_file","quiet");
	&run_cmd("rm -rf $tmp_dir/*.*");
}



sub run_signalp5{
	my $in_file=$_[0];
	my $model=$_[1];
	my $out_file=$_[2];
	my $final_fasta=$_[3];
	my $cwd = getcwd();
	$in_file=$cwd."/".$in_file;
	my $prefix=$cwd."/Tmp";
	# $final_fasta=$cwd."/".$final_fasta;
	# $out_file=$cwd."/".$out_file;
	chdir($path_signalp);
	&run_cmd("./signalp -format short -mature -org $model -fasta $in_file -prefix $prefix");
	chdir($cwd)
	&run_cmd("mv Tmp_summary.signalp5 $out_file");
	&run_cmd("mv Tmp_mature.fasta $final_fasta");
	# &run_cmd("mv Tmp_mature.fasta $final_fasta");
}

sub run_signalp{
	my $in_file=$_[0];
	my $model=$_[1];
	my $out_file=$_[2];
	my $final_fasta=$_[3];
	&run_cmd("echo '' > $out_file");
	&run_cmd("echo '' > $final_fasta");
	my $tmp_fasta=$wdir."Tmp.fasta";
	my $out_tmp=$wdir."out_signalp_temp";
	my $out_tmp_fasta=$wdir."out_signalp_temp_matured.fasta";
	my $batch=1000;
	my $s1;
	open $s1,">",$tmp_fasta;
	open my $fa,"<",$in_file;
	my $c_c="";
	my $n=0;
	while(<$fa>){
		chomp($_);
		if ($_=~/^>/){
			$n++;
			if ($n==$batch){
				close $s1;
 				&run_cmd("$path_signalp -t $model -m $out_tmp_fasta $tmp_fasta > $out_tmp");
				&run_cmd("cat $out_tmp >> $out_file");
				if (-e $out_tmp_fasta){&run_cmd("cat $out_tmp_fasta >> $final_fasta");}
				&run_cmd("rm $out_tmp_fasta $out_tmp");
				open $s1,">",$tmp_fasta;
				$n=1;
			}
			print $s1 "$_\n";
		}
		else{
			print $s1 "$_\n";
		}
	}
	close $fa;
	close $s1;
	&run_cmd("$path_signalp -t $model -m $out_tmp_fasta $tmp_fasta > $out_tmp");
	&run_cmd("cat $out_tmp >> $out_file");
	if (-e $out_tmp_fasta){&run_cmd("cat $out_tmp_fasta >> $final_fasta");}
	&run_cmd("rm $out_tmp_fasta $out_tmp");
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
