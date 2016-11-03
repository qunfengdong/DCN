Electronic Medical Record (EMR) Network Analyser
------------------------------------------------
EMR Network Analyser is a pipeline of disease correlation analysis with retrospective matched cohort study design using Cox Proportional Hazards (Cox-PH) regression in combination of interactive network display using graph theory. It allows clinicians to explore the relationships between any statistically disease pairs easily by studying the network with customized filtering, rearranging the network and calculating all the possible path between disease pairs. 

## Prerequisities
* R >= 3.3.0

## Citation

In press. Will update when accepted.

## Install
To check out the source code, go to https://github.com/qunfengdong/EMR_network.git. To obtain the scripts and example fasta files, do the following:

```shell
$ git clone https://github.com/qunfengdong/EMR_network.git
```

After the github repository is cloned, you will find a folder named EMR_network. All the scripts and example data files will be included in it. 

## Quick start

This suite of analysis include four major parts, 1). Find all the exposed disease pairs; 2). Find all matched non-exposed disease pairs; 3). Perform Cox-PH regression on the two cohorts; 4). Visualize the results.

### Input file
* Before running any analysis using this pipeline, please make sure you have an input EMR file formatted as tab delimited with 6 columns in the order of patient ID, disease ID/name, disease diagnose date ([%Y-%m-%d] or [%Y/%m/%d]), age, gender, race. No header is needed. It should be something like the following:

```
67      17      2012-10-15      74      F       White
128     17      2014-05-11      24      F       White
138     9       2008-08-23      7       M       White
144     7       2008-05-22      73      F       White
144     8       2008-07-11      73      F       White
144     1       2008-07-25      73      F       White
144     25      2011-05-05      76      F       White
148     10      2008-11-30      97      F       White
148     19      2008-11-30      97      F       White
148     20      2008-11-30      97      F       White
148     2       2009-01-19      97      F       White
```

* If you have a seperate file containing the acutual name of each disease, you could supply it in the scripts, and the output would show the disease name instead of disease ID.

## Getting started

### Step 1
* Find all the exposed disease pairs. Intermediate files with \<Disease A\>.vs.\<Disease B\></Dis>.csv will be outputted for each pair.
```
$ Rscript a1_exposure.r -i test.tsv
```
More options available:
```
$ Rscript a1_exposure.r -h
doParallel    foreach   optparse 
      TRUE       TRUE       TRUE 
Usage: Rscript a1_exposure.r -i <filename> [options]

This R script will find the exposure population for all diesease pairs appearing in the input file that satisfy disease A -> disease B, and outputs CSV files of selected exposed patients for each pair of disease combinations. The output from this step will be the input for next step.

Options:
	-i CHARACTER, --infile=CHARACTER
		Input file name [required]. 
		Should be tab delimited with 6 columns in the order of patient_ID, Disease, Disease_Date ([%Y-%m-%d] or [%Y/%m/%d]), Age, Gender, Race. No header is needed. You could provide more information about the patient (additional columns), but only the measurement at disease A will be recorded

	-p INTEGER, --processors=INTEGER
		Number of cores/CPUs to use, default is 8

	-h, --help
		Show this help message and exit
```
During the process of setting up the database, NCBI's 16SMicrobial.tar.gz file, taxdmp.zip, and taxdb.tar.gz will be downloaded into a default folder: ./db/, and uncompressed. 16SMicrobial.taxID.taxonomy under the ./db directory is the taxonomy file should be supplied to the 2.blca.py as the database. And an environmental variable called BLASTDB has to be set up manually. There will be instruction at the end of this script to let you know what shell command you should run to set it up. **If you login to a new shell, this environmental variable has to be set up again before running the analysis.**

```
export BLASTDB=/location/of/taxdb.bti/and/taxdb.btd/
```

Normally, it should be located in the ./db/ directory.

### Step 2 
* Run your analysis with the compiled database. Please run:
```
$ python 2.blca.py -i test.fasta
```
If you are running your analysis somewhere else other than in the BLCA_script directory, please do the following:
```
$ python /location/to/2.blca.py -i test.fasta -r /location/to/your/database/16SMicrobial.taxID.taxonomy
```
More options are the following:
```
$ python 2.blca.py -h

<< Bayesian-based LCA taxonomic classification method >>

   Please make sure the following softwares are in your PATH:
	1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.
	2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
	3.Biopython should be installed locally.

Usage: python 2.blca.py -i <fasta file> [option]


Arguments:
 - Required:
	-i		Input fasta file.
 - Taxonomy Profiling Options [filtering of hits]:
	-n		Number of times to bootstrap. Default: 100
	-t		Extra number of nucleotides to include at the beginning and end of the hits. Default: 10
	-d		Proportion of hits to include from top hit. Default: 0.1 [0-1]
	-e		Minimum evalue to include for blastn. Default: 0.1
	-a		Minimum bitscore to include for blastn hits. Default: 100
	-c		Minimum coverage to include. Default: 0.95 [0-1]
	-b		Minimum identity score to include. Default: 95 [0-100]
	-r		Reference Taxonomy file for the Database. Default: db/16SMicrobial.taxID.taxonomy
	-q		Refernece blast database. Default: db/16SMicrobial
	-o		Output file name. Default: <fasta>.blca.out
 - Alignment Options:
	-m		Alignment match score. Default: 1
	-f		Alignment mismatch penalty. Default: -2.5
	-g		Alignment gap penalty. Default: -2
 - Other:
	-h		Show program usage and quit
```

## Output
* A text file with sequence id in the first column, and taxonomy annotation with confidence scores after each level of annotaion (superkingdom, phylum, class, order, family, genus, species, strain).
* If no taxonomy annotation is available, it is listed as 'Not Available'

### Example output file:
```
Read1      superkingdom:Bacteria;100.0;phylum:Proteobacteria;100.0;class:Gammaproteobacteria;100.0;order:Enterobacterales;100.0;family:Enterobacteriaceae;100.0;genus:Lelliottia;100.0;species:Lelliottia nimipressuralis;100.0;strain:Lelliottia nimipressuralis;100.0;
Read2     superkingdom:Bacteria;100.0;phylum:Bacteroidetes;100.0;class:Cytophagia;100.0;order:Cytophagales;100.0;family:Cytophagaceae;100.0;genus:Runella;100.0;species:Runella slithyformis;100.0;strain:Runella slithyformis DSM 19594;97.5;
Read3      superkingdom:Bacteria;100.0;phylum:Actinobacteria;100.0;class:Actinobacteria;100.0;order:Streptomycetales;100.0;family:Streptomycetaceae;100.0;genus:Streptomyces;100.0;species:Streptomyces echinatus;100.0;strain:Streptomyces echinatus;50.0;
```

## Version
* Version 1.2 An alternative public release

## Authors
* Dr. Xiang Gao, theoretical conception and algorithm development
* Dr. Qunfeng Dong, algorithm development
* Huaiying Lin, program coding and testing
* Kashi Revanna, program coding and package development

## Error report

Please report any errors or bugs to hlin2@luc.edu.

## License
GNU

## Acknowledgements