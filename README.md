Electronic Medical Record (EMR) Network Analyser
------------------------------------------------
EMR Network Analyser is a pipeline of disease correlation analysis with retrospective matched cohort study design using Cox Proportional Hazards (Cox-PH) regression in combination of interactive network display using graph theory. It allows clinicians to explore the relationships between any statistically disease pairs easily by studying the network with customized filtering, rearranging the network and calculating all the possible path between disease pairs. 

## Prerequisities

* Linux
* python >= 2.7.0
* Rscript >= 3.3.0

## Citation

In press. Will update when accepted.

## Install

To check out the source code, go to https://github.com/qunfengdong/EMR_network. To obtain the scripts and example EMR files, do the following:

```shell
$ git clone https://github.com/qunfengdong/EMR_network.git
```

After the github repository is cloned, you will find a folder named EMR_network. All the scripts and example data files will be included in it. 

## Quick start

This suite of analysis include four major parts:

1. Find all the exposed disease pairs; 
2. Find all matched non-exposed disease pairs; 
3. Perform Cox-PH regression on the two cohorts; 
4. Visualize the results;
5. Optional: Validate a certain disease pair using Random Forest survival analysis.

Step 1-3 is achieved by a1_runEMR.r, step 4 is achieved by a2_parse.py and index.html, step 5 is achieved by a3_RFopt.r.

### Required input files

* Before running any analysis using this pipeline, please make sure you have an input EMR file formatted as tab delimited with 6 columns in the order of patient ID, disease ID/name, disease diagnose date ([%Y-%m-%d] or [%Y/%m/%d]), age, gender, race. No header is needed. It should be something like the following:

```
P00001  62      2012-01-16      66      F       White
P00002  509     2009-01-12      68      M       African American
P00002  23      2014-11-07      73      M       African American
P00003  23      2015-07-20      78      M       White
P00007  44      2014-04-28      51      M       White
P00008  509     2010-07-12      54      M       White
P00015  509     2008-04-21      60      M       White
P00016  63      2008-03-16      83      F       White
P00017  23      2014-11-13      48      F       White
P00018  509     2011-07-03      76      M       African American
P00018  14      2012-09-20      77      M       African American
P00018  23      2014-01-09      79      M       African American
P00019  23      2015-01-09      89      M       White
P00020  14      2012-05-25      65      F       African American
```

* Meta file should be tab separated with header. The first column is the disease ID (should be consistant with the disease ID you have in the input data file), and second column is the disease name. Example of meta file:

| id_col | name_long |
|--------|-----------|
| 14     | D014      |
| 23     | D023      |
| 29     | D029      |
| 509    | D509      |
| 43     | D043      |
| 44     | D044      |
| 59     | D059      |
| 62     | D062      |
| 63     | D063      |
| 67     | D067      |

### Test dataset

* A test dataset with those two files are availabe in test.tar.gz file. Please unzip it to view and go through the tutorial.

```
$ tar zxvf test.tar.gz
```

1. View the example input file: cdc_first_test.tsv
2. View the example meta file: disease_code_test.tsv
3. View the intermediate time table: \<DiseaseA\>.vs.\<DiseasB\>.csv
4. View the survival curves of significant disease pairs based on Cox-PH regression: \<DiseaseA\>_\<DiseasB\>.png
5. Click on the index.html to view the result from example data.
6. View the survival curves of significant disease pairs based on Random Forest survival analysis: \<DiseaseA\>.vs.\<DiseasB\>.RFsurv.png

* Note: If you re-run the tutorial, all present files will be overwritten. If you want to compare your results to the provided example results, please save the files 3-6 in a newly created folder.

## Getting started

### Step 1: a1_runEMR.r 
* Find all the exposed and matching non-exposed disease pairs. And Cox-PH regression will be performed taking age, race, gender as covariates (or any other confounding factors supplied by user), and final output is an adjacency matrix with all the adjusted significant hazard ratios and survival curve graphs. Intermediate files containing the exposed and non-exposed cohorts with name as \<Disease A\>.vs.\<Disease B\></Dis>.csv can be outputted for each disease pair if requested with flag "-r". Those intermediate files are required for a3_RFopt.r step as the input file. **Turn on the flag if you want to validate the results from a1_runEMR.r**. 

```
$ Rscript a1_runEMR.r -i cdc_first_test.tsv -m disease_code_test.tsv
```

More options available:

```
$ Rscript a1_runEMR.r -h
doParallel    foreach     OIsurv   optparse    ggplot2  ggfortify 
      TRUE       TRUE       TRUE       TRUE       TRUE       TRUE 
Usage: Rscript a1_runEMR.r -i <inputFileName> -m <metaFileName> [options]

This R script will first find the exposure population for all diesease pairs appearing in the input file that satisfy disease A -> disease B, and and second, it will find the matching non-exposed population for all diesease pairs with some criteria, such as the matching patient should have the same gender and race as control. After the non-exposed are matched with the exposed group, Cox-PH regression will be performed taking age, race, gender as covariates (or any other confounding factors supplied by user), and final output is an adjacency matrix with all the significant hazard ratios, and survival curve graphs. Multiple test correction will be applied to the p-values.

Options:
	-i CHARACTER, --infile=CHARACTER
		Input file name [required]. 
		Should be tab delimited with 6 columns in the order of patient_ID, Disease, Disease_Date ([%Y-%m-%d] or [%Y/%m/%d]), Age, Gender, Race. No header is needed. You could provide more information about the patient (additional columns), but only the measurement at disease A will be recorded

	-m CHARACTER, --metafile=CHARACTER
		Meta file for each disease name [required]. 
		First column is the disease ID used in the input file column Disease, second column is the disease name. All the other columns should be additional attribute of the disease, tab separated. Header is needed

	-o CHARACTER, --outfile=CHARACTER
		Adjacency matrix output file name, default is all.edges.csv

	--outfolder=CHARACTER
		Intermediate files and survival PNG output folder, default is current directory

	--maxn=INTEGER
		Maximum number of subjects to be qualified as valid exposed cohort, default is 5000. If the total number of patients exceed this number, a subset of 5000 subjects will be randomly selected.

	--minn=INTEGER
		Minimum number of subjects to be qualified as valid exposed cohort, default is 100

	-c INTEGER, --caseNumber=INTEGER
		Minimum number of A->B cases to keep disease pair, default is 20

	--duration=INTEGER
		study duration (years), default is 10 years

	-y INTEGER, --year=INTEGER
		Age difference between exposed and non-exposed, default is 5 years

	-d INTEGER, --days=INTEGER
		Visit date difference between exposed and non-exposed, unit: days, default is 7 days

	-s DOUBLE, --sigcut=DOUBLE
		Significant P-value cutoff. Only p-values that are less than this cutoff will be considered as significant P-values, and kept in the adjacency matrix. Default is 0.05

	-e DOUBLE, --exp_coef=DOUBLE
		Exponentiated coefficients cutoff used in plotting survival curve. Default is 1

	-a CHARACTER, --adjust=CHARACTER
		Method of adjusting p-values, choices are holm, hochberg, hommel, bonferroni, BY, and fdr. The default is BH

	-r, --writeintermediates
		Write out individual disease pair survival tables, default is FALSE

	-p INTEGER, --processors=INTEGER
		Number of cores/CPUs to use, default is 8

	-h, --help
		Show this help message and exit
```

* Output: 
1. all.edges.csv \[default name\]: adjacency matrix.
2. Survival curve graphs if tested significant.
3. (Optional) Time table for each disease pair.

* Example **all.edges.csv** file (This won't match the result after you run the test.tsv):


| from | to  | coef        | exp_coef    | se_coef     | z           | Pr          | N    | HRtest.p    | Pr_adjusted | from_name_long | to_name_long |
|------|-----|-------------|-------------|-------------|-------------|-------------|------|-------------|-------------|----------------|--------------|
| 14   | 23  | 1.427457704 | 4.168089196 | 0.095763325 | 14.90609999 | 0           | 2746 | 0.999988358 | 0           | D014           | D023         |
| 509  | 23  | 0.683504686 | 1.980807691 | 0.084860393 | 8.054460489 | 7.77E-16    | 2896 | 0.552183633 | 8.16E-15    | D509           | D023         |
| 63   | 23  | 0.662536077 | 1.939705344 | 0.145536957 | 4.552356253 | 5.30E-06    | 830  | 0.991949308 | 2.48E-05    | D063           | D023         |
| 63   | 62  | 2.911744586 | 18.38885154 | 0.600664408 | 4.847539734 | 1.25E-06    | 830  | 0.999999995 | 6.56E-06    | D063           | D062         |
| 43   | 23  | 0.855137206 | 2.351697024 | 0.130494105 | 6.553071524 | 5.64E-11    | 1316 | 0.422697221 | 4.73E-10    | D043           | D023         |
| 29   | 23  | 0.842790839 | 2.322840614 | 0.16442774  | 5.125600098 | 2.97E-07    | 594  | 0.754490443 | 1.78E-06    | D029           | D023         |
| 23   | 509 | 1.486126685 | 4.419942492 | 0.114816438 | 12.94350102 | 0           | 6388 | 0.243396393 | 0           | D023           | D509         |
| 23   | 14  | 1.011870161 | 2.750740535 | 0.097769744 | 10.34952241 | 0           | 6388 | 0.932861154 | 0           | D023           | D014         |
| 23   | 63  | 0.481728007 | 1.618869405 | 0.164242912 | 2.933021586 | 0.003356806 | 6388 | 0.990998063 | 0.010070417 | D023           | D063         |
| 23   | 43  | 0.739185893 | 2.094229893 | 0.136233051 | 5.425892526 | 5.77E-08    | 6388 | 0.99847771  | 4.04E-07    | D023           | D043         |
| 23   | 44  | 0.574941042 | 1.777025754 | 0.189093351 | 3.040514325 | 0.002361745 | 6388 | 0.848554064 | 0.007630252 | D023           | D044         |
| 23   | 29  | 1.045118223 | 2.8437347   | 0.296102045 | 3.529587994 | 0.000416207 | 6388 | 0.999983025 | 0.001589155 | D023           | D029         |
| 62   | 63  | 1.624442928 | 5.075590779 | 0.49812411  | 3.26112086  | 0.001109727 | 338  | 0.998567186 | 0.003884046 | D062           | D063         |

* Example survial curve: 

<a href="https://github.com/yingeddi2008/scratch/blob/master/14_23.png"><img src="https://github.com/yingeddi2008/scratch/blob/master/14_23.png" alt="example of survival curve" style="float: right" width="450"/></a>

* Example time table:

| ID     | type | time | event | V3       | V4 | V5 | V6               |
|--------|------|------|-------|----------|----|----|------------------|
| P41812 | 0    | 350  | 0     | 9/4/11   | 68 | F  | White            |
| P39443 | 0    | 770  | 0     | 11/7/10  | 52 | M  | White            |
| P41112 | 0    | 0    | 0     | 3/4/10   | 70 | M  | White            |
| P37186 | 0    | 0    | 0     | 12/24/16 | 67 | F  | White            |
| P28958 | 0    | 146  | 0     | 10/5/14  | 59 | F  | White            |
| P15618 | 0    | 1004 | 0     | 6/12/08  | 67 | M  | White            |
| P14117 | 0    | 0    | 0     | 9/12/08  | 78 | F  | White            |
| P24383 | 0    | 125  | 0     | 9/1/08   | 20 | F  | White            |
| P00008 | 1    | 0    | 0     | 9/4/15   | 60 | M  | White            |
| P00010 | 1    | 0    | 0     | 9/18/11  | 79 | M  | White            |
| P00018 | 1    | 0    | 0     | 1/9/14   | 79 | M  | African American |
| P00028 | 1    | 805  | 0     | 6/5/10   | 55 | M  | White            |
| P00040 | 1    | 222  | 0     | 1/4/10   | 85 | F  | White            |

### Step 2: a2_parse.py
* Parse the adjacency matrix into a json object for webpage network display.

```
$ python a2_parse.py <output file: all.edges.csv from a1_runEMR.r>
```

* Then just click-open the index.html, you will see the network as shown similar to the following.

* A screenshot of final network display.

<a href="https://github.com/yingeddi2008/scratch/blob/master/EMR_display.png"><img src="https://github.com/yingeddi2008/scratch/blob/master/EMR_network.PNG" alt="Screenshot of network disply"/></a>

### Step 3
* Visualization with cytoscape.js. The output from step 2 needs to be formatted into a json object before viewing.

```
$ python a3_parse.py
```


## Version
* Version 0.9 

## Authors

* Huaiying Lin. M.S., algorithm development, program coding and testing
* Dr. Xiang Gao, theoretical conception and algorithm development
* Dr. Qunfeng Dong, algorithm development
* Kashi Revanna, cytoscape visualization

## Error report

Please report any errors or bugs to hlin2@luc.edu.

## License
GNU
