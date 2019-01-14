# Inovirus_detector
This set of script can be used to identify putative inovirus sequences in draft genome assemblies or metagenome assemblies. 

## Requirements
* Bioperl, blast, hmmer, and randomforest from R, all available through conda
```
conda create -n inovirus_detector -c bioconda perl-bioperl=1.7.2 blast=2.7.1 hmmer=3.2.1 
conda install -n inovirus_detector -c r r-randomforest
```
* SignalP and TMHMM, not available through conda.
Download and install SignalP according to instructions from the SignalP v4.1 authors (http://www.cbs.dtu.dk/services/doc/signalp-4.1.readme)
Download and install Tmhmm according to instructions from the Tmhmm v2.0 authors (http://www.cbs.dtu.dk/services/doc/tmhmm-2.0c.readme)

* Extract the Inovirus_db.tar.gz archive
```
tar -xvf Inovirus_dbs.tar.gz
```

* PFAM database for de novo annotation (note: if this database is already in your system, you can provide the path as an argument to Step 1):
```
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
```

## Step 1: De novo identification of putative inovirus sequences in a contig set - input: (meta)genome contig(s) in genbank format
### Input example: 
Example_files/2731957639/2731957639_129103.assembled.gbk
```
source activate inovirus_detector
./Identify_candidate_fragments_from_gbk.pl -g Example_files/2731957639/2731957639_129103.assembled.gbk -p Pfam-A.hmm
```
### This first step will
* Look for pI-like proteins in the genome provided in the genbank file
* Perform a custom and simple annotation for 30kb windows around the putative pI-like protein(s)
* Generates for each putative inovirus fragment (i.e. 30kb around a putative pI-like protein) a set of files including an annotation gff file and the corresponding nucleotide fragment in fasta file (these are the two input files for Step 2).
Candidate fragments and their associated files are identified using the name of the protein initially detected as a putative pI-like protein.
For example: when run on example file CP007542_Synechocystis_PCC_6714.gb, the result files should include a gff of candidate inovirus fragments similar to Example_expected_results/CP007542_Synechocystis_PCC_6714_frag_D082_13100_annot.gff and a fasta file similar to Example_expected_results/CP007542_Synechocystis_PCC_6714_frag_D082_13100_nucl.fna

### Note:
If error: "ListUtil.c: loadable library and perl binaries are mismatched", this is a known conda issue, that can be fixed with the following steps:
Create a file etc/conda/activate.d/update_perllib.sh in your conda environment folder including the following lines:
```
#!/bin/sh
export OLD_PERL5LIB=$PERL5LIB
export PERL5LIB=`pwd`/../../../lib/site_perl/5.26.2/
```
Then create a file etc/conda/deactivate.d/update_perllib.sh in your conda environment folder including the following lines:
```
#!/bin/sh
export PERL5LIB=$OLD_PERL5LIB
```
### Note:
This custom annotation uses tmhmm and signalp to identify putative inovirus coat proteins. The same can be done on individual protein fasta file using the stand-along script "Predict_inovirus_coat_proteins.pl"

## Step 2: Refine inovirus genome prediction - input: gff annotation of putative fragments, either from step 1 or from a custom annotation pipeline. 
### Example input files: 
Example_files/2731957639_expected_results/2732535622_129103.assembled_frag_Ga0128599_102362_annot.gff & Example_files/2731957639_expected_results/2732535622_129103.assembled_frag_Ga0128599_102362_nucl.fna (output from step 1)
For candidate fragment 1
```
./Get_inovirus_prediction_score_from_gff_fragments.pl -i Example_files/2731957639_expected_results/2732535622_129103.assembled_frag_Ga0128599_102362_annot.gff -f Example_files/2731957639_expected_results/2732535622_129103.assembled_frag_Ga0128599_102362_nucl.fna
```
For candidate fragment 2
```
./Get_inovirus_prediction_score_from_gff_fragments.pl -i Example_files/2731957639_expected_results/2732535622_129103.assembled_frag_Ga0128599_10338_annot.gff -f Example_files/2731957639_expected_results/2732535622_129103.assembled_frag_Ga0128599_10338_nucl.fna
```
Final result files should be identical to 2732535622_129103.assembled_frag_Ga0128599_102362_annot_inovirus-predictions.csv, 2732535622_129103.assembled_frag_Ga0128599_102362_annot_inovirus-predictions-refined.csv, 2732535622_129103.assembled_frag_Ga0128599_10338_annot_inovirus-predictions.csv, 
and 2732535622_129103.assembled_frag_Ga0128599_10338_annot_inovirus-predictions-refined.csv in folder Example_files/2731957639_expected_results/

### This second step will
* Run the random forest classifier to identify which subset of the candidate fragment if most likely an inovirus genome (if any)
* Attempt to refine the prediction by detecting canonical attachment (att) site, i.e. direct repeats flanking a provirus with one repeat in a tRNA or outside of an integrase
* Otherwise, attempt to refine the prediction by detecting non-canonical attachment (att) site, i.e. direct repeats neither in a tRNA nor outside of an integrase

### Note
If using a gff from a custom annotation pipeline (i.e. not the script used in Step 1), formatting should follow the following rules:
Each fragment should be preceeded by a line starting with \"## Fragment_id\", The Fragment_id must be the gene id of the putative morphogenesis gene
PFAM annotation should be indicated as \"pfam=PFAM_domain_name Score\" in the 9th field
Inoviridae PC annotation with score should be indicated as \"inoPC=PC_XXX Score\" in the 9th field as well
Prediction of coat protein should be indicated as \"Coat_pred=yes/no\" in the 9th field as well
Note: the fragment is expected to be a microbial genome fragment, and start with a CDS. 
The starting coordinate of the first CDS will thus be used to \"shift\" the coordinate when dealing with the fragment sequence (i.e. transform the genes coordinates from the global genome reference to the fragment reference)
if the starting coordinate of the first CDS should not be used as the start of the fragment, you can indicate in the header line after the fragment id (tab separated) a number of bp to use in this transformation (can be 0 if the fragment is an entire contig)
