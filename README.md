Electronic Medical Record (EMR) Network Analyser
------------------------------------------------
EMR Network Analyser is a pipeline of disease correlation analysis with retrospective matched cohort study design using Cox Proportional Hazards (Cox-PH) regression in combination of interactive network display using graph theory. It allows clinicians to explore the relationships between any statistically disease pairs easily by studying the network with customized filtering, rearranging the network and calculating all the possible path between disease pairs. 

## Prerequisities

* Linux
* R >= 3.3.0

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
4. Visualize the results.

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
* Find all the exposed disease pairs. Intermediate files with \<Disease A\>.vs.\<Disease B\></Dis>.csv will be outputted for each pair. This step completes the first part of this pipeline.
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

### Step 2 
* According to the exposed population identified from the first step, matched non-exposed population will be chosen. After the non-exposed are matched with the exposed, Cox-PH regression will be performed taking age, race, gender as covariates (or any other confounding factors supplied by user), and final output is an adjacency matrix with all the adjusted significant hazard ratios and survival curve graphs.

```
$ Rscript a2_non_exposed.cox.adjust.r -i test.tsv
```

More options are the following:
```
$ Rscript a2_non_exposed.cox.adjust.r -h
doParallel    foreach     OIsurv   optparse 
      TRUE       TRUE       TRUE       TRUE 
Usage: Rscript a2_non_exposed.cox.adjust.r -i <filename> [options]

This R script will find the matching non-exposed population for all diesease pairs appearing in the output file from last step with some criteria, such as the matching patient should have the same gender and race as the exposed subject. After the non-exposed are matched with the exposed, Cox-PH regression will be performed taking age, race, gender as covariates (or any other confounding factors supplied by user), and final output is an adjacency matrix with all the adjusted significant hazard ratios, and survival curve graphs.

Options:
	-i CHARACTER, --infile=CHARACTER
		Input file name [required]
		Should be tab delimited with 6 columns in the order of patient_ID, Disease, Disease_Date ([%Y-%m-%d] or [%Y/%m/%d]), Age, Gender, Race. No header is needed. You could provide more information about the patient (additional columns), but only the measurement at disease A will be recorded. And these/this additional variable(s) will be used as covariate during Cox-PH regression, but won't be used as matching criteria

	-o CHARACTER, --outfile=CHARACTER
		Adjacency matrix output file name, default is all.edges.csv

	-m CHARACTER, --metafile=CHARACTER
		Meta file for each disease name. First column is the disease ID used in the input file column Disease, all the other column should be additional attribute of the disease, tab separated. Header is needed

	-a CHARACTER, --adjust=CHARACTER
		Method of adjusting p-values, choices are holm, hochberg, hommel, bonferroni, BY, and fdr. The default is BH

	-s DOUBLE, --sigcut=DOUBLE
		Significant P-value cutoff. Only p-values that are less than this cutoff will be considered as significant P-values, and kept in the adjacency matrix. Default is 0.05

	-y INTEGER, --year=INTEGER
		Age difference between exposed and non-exposed, default is 5 years

	-d INTEGER, --days=INTEGER
		Visit date difference between exposed and non-exposed, unit days, default is 7 days

	-p INTEGER, --processors=INTEGER
		Number of cores/CPUs to use, default is 8

	-h, --help
		Show this help message and exit
```

### Step 3
* Visualization with cytoscape.js.

## Output
* A CSV file with sequence id in the first column, and taxonomy annotation with confidence scores after each level of annotaion (superkingdom, phylum, class, order, family, genus, species, strain).
* If no taxonomy annotation is available, it is listed as 'Not Available'

### Example output file:

| from | to | coef               | exp_coef         | se_coef            | z                | Pr                   | N     | Pr_adjusted          | from_cdc_desc | to_cdc_desc                                  | 
|------|----|--------------------|------------------|--------------------|------------------|----------------------|-------|----------------------|---------------|----------------------------------------------| 
| 10   | 18 | 0.772471921293428  | 2.16511162992917 | 0.112934793120933  | 6.83998172703295 | 7.92033105767587e-12 | 16064 | 8.59355919757832e-11 | Depression    | Alzheimer Disease                            | 
| 10   | 11 | 0.277057388828638  | 1.31924207862013 | 0.0338307755002824 | 8.18950747452789 | 2.22044604925031e-16 | 15952 | 3.21224528458212e-15 | Depression    | Rheumatoid Arthritis/Osteoarthritis          | 
| 10   | 19 | 0.706414526729088  | 2.02671149585188 | 0.0749407576526056 | 9.42630617645652 | 0                    | 16250 | 0                    | Depression    | Senile Dementia                              | 
| 10   | 20 | 0.0835078464368838 | 1.08709374507187 | 0.0313596169960636 | 2.66291027876285 | 0.00774680666244776  | 16032 | 0.0235112873531631   | Depression    | Anemia                                       | 
| 10   | 31 | 0.111767501241316  | 1.11825283802201 | 0.0476239027973893 | 2.34687824130707 | 0.0189314401604171   | 16010 | 0.0490846030097312   | Depression    | Liver Disease                                | 
| 10   | 28 | 0.266936020103298  | 1.30595688867084 | 0.0387219161535401 | 6.8936676337205  | 5.43720624079924e-12 | 16012 | 6.0506346371971e-11  | Depression    | Urinary Incontinence                         | 
| 10   | 33 | 0.510444628015953  | 1.66603179469938 | 0.172117419224065  | 2.96567674740376 | 0.0030201774719375   | 16184 | 0.0100057788001594   | Depression    | Parkinson's Disease                          | 
| 10   | 32 | 0.132298841532258  | 1.14144938078497 | 0.0504671483452885 | 2.62148438875701 | 0.00875477718516648  | 15988 | 0.0261139058306684   | Depression    | Sleep Apnea                                  | 
| 10   | 5  | 0.210064999380328  | 1.23375825087232 | 0.0600064275096725 | 3.50070830906352 | 0.000464023431176241 | 16120 | 0.00179809079580793  | Depression    | Chronic Obstructive Pulmonary Disease (COPD) | 
| 10   | 9  | 0.241120268657163  | 1.27267408890795 | 0.0521925903589008 | 4.61981800479927 | 3.84076791992705e-06 | 16208 | 2.15083003515915e-05 | Depression    | Asthma                                       | 

## Version
* Version 0.9 An alternative public release

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