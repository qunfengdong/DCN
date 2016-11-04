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
* Find all the exposed disease pairs. Intermediate files with \<Disease A\>.vs.\<Disease B\></Dis>.csv will be outputted for each disease pair. Those intermediate files are required for next step. **Do not remove them**. This step completes the first part of this pipeline.
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
* According to the exposed population identified from the first step (In each of the CSV file), matched non-exposed population will be chosen. After the non-exposed are matched with the exposed, Cox-PH regression will be performed taking age, race, gender as covariates (or any other confounding factors supplied by user), and final output is an adjacency matrix with all the adjusted significant hazard ratios and survival curve graphs. This script completes the second and third parts of this pipeline.

```
$ Rscript a2_non_exposed.cox.adjust.r -i test.tsv -m test.meta.tsv
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

* Example of meta file:

| ID | Name          | 
|----|---------------| 
| 1  | Heart Failure | 
| 3  | Liver cancer  | 
| 4  | Diabetes      | 
| 5  | COPD          | 
| 7  | Hypertension  | 
| 9  | Allergy       | 
| 10 | Depression    | 

* The output of this step are the following: 

1. Intermediate files named as &lt;Disease A&gt;.vs.&lt;Disease B&gt;.all.csv. You can look into the raw input file for a specific disease pair of your interest.
2. Survival curves PNG files named as &lt;Disease A&gt;.vs.&lt;Disease B&gt;.png.
3. CSV file named **all.edges.csv** with the _from_ disease ID, _to_ disease ID, coefficient, exponentiated coefficient, standard error of the coefficient, z value, p-value, sample size (N), adjusted p-value, the _from_ disease name and _to_ disease name as columns. **If you do not supply a meta file when running, the _from_ disease name and _to_ disease name won't appear in the all.edges.csv.**

* Example survial curve:

![alt text][display]
[display]: https://github.com/yingeddi2008/scratch/blob/master/4_7.png "Example survival curve graph"

* Example **all.edges.csv** file (This won't match the result after you run the test.tsv):


| from | to | coef  | exp_coef | se_coef | z     | Pr    | N     | Pr_adjusted | from_Name | to_Name                                  | 
|------|----|-------|----------|---------|-------|-------|-------|-------------|---------------|----------------------------------------------| 
| 10   | 18 | 0.772 | 2.165    | 0.113   | 6.840 | 0.000 | 16064 | 8.6E-11     | Depression    | Alzheimer Disease                            | 
| 10   | 11 | 0.277 | 1.319    | 0.034   | 8.190 | 0.000 | 15952 | 3.2E-15     | Depression    | Rheumatoid Arthritis/Osteoarthritis          | 
| 10   | 19 | 0.706 | 2.027    | 0.075   | 9.426 | 0.000 | 16250 | 0.0E+00     | Depression    | Senile Dementia                              | 
| 10   | 20 | 0.084 | 1.087    | 0.031   | 2.663 | 0.008 | 16032 | 2.4E-02     | Depression    | Anemia                                       | 
| 10   | 31 | 0.112 | 1.118    | 0.048   | 2.347 | 0.019 | 16010 | 4.9E-02     | Depression    | Liver Disease                                | 
| 10   | 28 | 0.267 | 1.306    | 0.039   | 6.894 | 0.000 | 16012 | 6.1E-11     | Depression    | Urinary Incontinence                         | 
| 10   | 33 | 0.510 | 1.666    | 0.172   | 2.966 | 0.003 | 16184 | 1.0E-02     | Depression    | Parkinson's Disease                          | 
| 10   | 32 | 0.132 | 1.141    | 0.050   | 2.621 | 0.009 | 15988 | 2.6E-02     | Depression    | Sleep Apnea                                  | 
| 10   | 5  | 0.210 | 1.234    | 0.060   | 3.501 | 0.000 | 16120 | 1.8E-03     | Depression    | Chronic Obstructive Pulmonary Disease (COPD) | 
| 10   | 9  | 0.241 | 1.273    | 0.052   | 4.620 | 0.000 | 16208 | 2.2E-05     | Depression    | Asthma                                       | 


### Step 3
* Visualization with cytoscape.js.

* A screen shot of final network display.

![alt text][display]
[display]: https://github.com/yingeddi2008/scratch/blob/master/EMR_display.png "Example Network Display"

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
