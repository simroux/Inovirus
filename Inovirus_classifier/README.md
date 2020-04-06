# Inovirus_detector
This script can be used to do a blast-based affiliation of new inovirus genomes

## Requirements
* blast, hmmer, and prodigal, all available through conda
```
conda create -n inovirus_classifier -c bioconda blast=2.9.0 prodigal=2.6.3 hmmer=3.3
```

## Using the tool
* After loading the conda environment, run ./Run_Ino_classifier.pl to get a list of parameters. Required parameters are:
- -i: path to a fasta file of inovirus genomes to be classified
- -d: path to the "Ino_classifier_db/" database
- -o: path to the desired output directory (will be created if needed)

## Example data
* The folder "Example_data/" has example input and output files including 3 genomes: 1 Inoviridae, 1 Plectroviridae, and 1 Paulinoviridae.
* Run as e.g.
```
./Run_Ino_classifier.pl -i Example_Data/Input/Input_genomes.fna -d Ino_classifier_db/ -o Example_Data/Output/
```
The content of the Example_Data/Output/ should be identical to the files in Example_Data/Expected_Output/
