Disease Correlation Network (DCN) Analyser
------------------------------------------------
Disease Correlation Network (DCN) Analyser is a pipeline of disease correlation analysis with retrospective matched cohort study design using Cox Proportional Hazards (Cox-PH) regression in combination of interactive network display using graph theory. It allows clinicians to explore the relationships between any statistically disease pairs easily by studying the network with customized filtering, rearranging the network and calculating all the possible path between disease pairs. 

## Prerequisities

* Linux or Mac OS
* python >= 3.7.0
* Rscript >= 3.5.0

## Environment setup

If you do not have the above environment, please set it up via the following <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/" target="_blank">conda</a> commands:

```bash
$ conda create --name DCN python=3.7 r=3.5 r-essentials
$ conda activate DCN
```
You can change the environment name `DCN` to anything you want.

After you are done, use the following code to deactivate.
```bash
$ conda deactivate
```

## Citation

In press. Will update when accepted.

Live interactive network display for the paper: <a href="http://cbi.lumc.edu/disease/" target="_blank">http://cbi.lumc.edu/disease/</a>.

## Install

To check out the source code, go to https://github.com/qunfengdong/DCN. To obtain the scripts and example DCN files, do the following:

```
$ git clone https://github.com/qunfengdong/DCN.git
$ cd DCN
```

After the github repository is cloned, you will find a folder named DCN. All the scripts and example data files will be included in it. 

### Test dataset

* A test dataset with those two files are available in test.tar.gz file. Please unzip it to view and go through the tutorial. _Please note that the test dataset is simulated, thus its results are not a true reflection of scientific discovery_.

```bash
$ tar zxvf test.tar.gz
```

1. View the example input file: **cdc_first_test.tsv**
2. View the example meta file: **disease_code_test.tsv**
3. View the time table: \<DiseaseA\>_\<DiseaseB\>.csv
4. View the survival curves of significant disease pairs based on Cox-PH regression: \<DiseaseA\>\_\<DiseaseB\>.cox.survival.png and residual plot: \<DiseaseA\>\_\<DiseaseB\>.cox.residual.png
5. Click on the index.html to view the result from example data.
6. View the survival curves of significant disease pairs based on Random Forest survival analysis: \<DiseaseA\>_\<DiseaseB\>.rf.survival.png

#### Running test

You will need to run the following codes inside the DCN folder to complete the test:

* Generate disease pairs
```bash
$ Rscript a0_generator.r -i ./test/cdc_first_test.tsv -m ./test/disease_code_test.tsv --outfolder ./example
```
* Analyze disease pairs using both Cox-PH and RF regression
```bash
$ Rscript a1_analyzer.r -m ./test/disease_code_test.tsv --inputfolder ./example --outputfolder ./example/output
$ Rscript a1_analyzer.r --method RF -m ./test/disease_code_test.tsv --inputfolder ./example --outputfolder ./example/output 
```
* Convert resulting Cox-PH regression results to json format for viewing
```bash
$ python a2_parse.py ./example/output/CoxPH.edge.csv
```
* Move results to the server folder
```bash
$ cp ./example/output/* web_server
```

Then just click-open the index.html inside web_server, you will see the network display. You can compare your results under the `example` folder with files in `test` folder.

**NOTE**: Files: `CoxPH.edge.csv`, `RF.edge.csv`, `all.edges.csv.js`, and all `PNG` files should be in the same folder as the `index.html`.

## Quick start

* This suite of analysis includes five major parts:

1. Find all the exposed disease pairs; 
2. Find all matched non-exposed disease pairs; 
3. Perform Cox-PH regression on the two cohorts; 
4. Perform Random Forest survival analysis on the two cohorts;
5. Visualize the results

* Step 1-2 is achieved by a0_generator.r, step 3-4 is achieved by a1_analyzer.r, and step 5 is achieved by a2_parse.py and index.html.

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

| id_col | name_long              |
|--------|------------------------|
| 14     | Prostate Cancer        |
| 23     | Glaucoma               |
| 29     | Obesity                |
| 509    | Med:Progestins         |
| 43     | Acute Myeloid Leukemia |
| 44     | Vitiligo               |
| 59     | Chalazion              |
| 62     | Cirrhosis            |
| 63     | Psoriasis              |
| 67     | Sarcoidosis            |

## Getting started

**NOTE**: You will need to enter your own file paths in the code below.

### Step 1: a0_generator.r
Find all the exposed and matching non-exposed disease pairs. 
* The default match finding criteria is: the matching subject is within 5 years of age difference (You can change this number using the argument: `-y`, `--year`), visit time is within 7 days (`-d`, `--days`) compared to the exposed counterparts, same race and gender. 
* The default screen criteria for a disease is: the minimum number of subjects for a disease is 100 (`--minn`), while the maximum is 10000 (`--maxn`), and the minimum number of disease A to disease B cases is 20 (`-c`, `--minCaseNumber`). 
* We have the option to choose whether to write out results for time-to-event = 0 (`--saveT0`). 

* Example code:
 
```bash
$ Rscript a0_generator.r -i /path/to/test/cdc_first_test.tsv -m /path/to/test/disease_code_test.tsv --outfolder /path/to/outputfolder
```

More options available:

```
$ Rscript a0_generator.r -h
data.table doParallel    foreach   optparse 
      TRUE       TRUE       TRUE       TRUE 
Usage: Rscript a0_generator.r -i <inputFileName> -m <metaFileName> [options]

		 >> This the step 1 of EMR package. << 
	This R script will find the exposure population for all diesease pairs appearing in the input file that satisfy disease A -> disease B, and and second, it will find the matching non-exposed population for all diesease pairs with some criteria, such as the matching patient should have the same gender and race as a good match control. Outputs are CSV files with paired exposed and non-exposed subjects information. One folder named T0_include with CSV files with time-to-event = 0 could be produced if you turn on the --saveT0 flag. Those files will be used as inputs for step 2 of the package.


Options:
	-i CHARACTER, --infile=CHARACTER
		Input file name [required]. 
		Should be tab delimited with 6 columns in the order of patient_ID, Disease, Disease_Date ([%Y-%m-%d] or [%Y/%m/%d]), Age, Gender, Race. No header is needed. You could provide more information about the patient (additional columns), but only the measurement at disease A will be recorded

	-m CHARACTER, --metafile=CHARACTER
		Meta file for each disease name [required]. 
		First column is the disease ID used in the input file column Disease, second column is the disease name. All the other columns should be additional attribute of the disease, tab separated. Header is needed

	--outfolder=CHARACTER
		Intermediate files output folder, default is current directory

	--maxn=INTEGER
		Maximum number of subjects to be qualified as valid exposed cohort, default is 10000. If the total number of patients exceed this number, a subset of 5000 subjects will be randomly selected.

	--minn=INTEGER
		Minimum number of subjects to be qualified as valid exposed cohort, default is 100

	-c INTEGER, --minCaseNumber=INTEGER
		Minimum number of A->B cases to keep disease pair, default is 20

	--duration=INTEGER
		study duration (years), default is 5 years

	-y INTEGER, --year=INTEGER
		Age difference between exposed and non-exposed, default is 5 years

	-d INTEGER, --days=INTEGER
		Visit date difference between exposed and non-exposed, unit: days, default is 7 days

	--saveT0
		A flag to save the cohort tables with time-to-event = 0. The output CSVs will be saved to T0_include folder, default is false

	-p INTEGER, --processors=INTEGER
		Number of cores/CPUs to use, default is 8

	-h, --help
		Show this help message and exit
```

* Output: \<DiseaseA\>_\<DiseaseB\>.csv files.

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


### Step 2: a1_analyzer.r
* Perform either Cox-PH regression or Random Forest survival analysis, taking age, race, gender as covariates (or any other confounding factors supplied by user), and final output is an adjacency matrix with all the adjusted significant P-values and survival curve graphs. 

#### Cox-PH regression
* The default is to perform Cox-PH regression. Example code:

```bash
$ Rscript a1_analyzer.r -m /path/to/disease_code_test.tsv --inputfolder /path/to/all_csv_files/ --outputfolder /path/to/all_csv_files/outputs
```

* Example **CoxPH.edge.csv** file output (This won't match the result after you run the code):

| from | to  | coef  | exp_coef | se_coef | z     | Pr    | N    | HRtest-p | Pr.BH | from_name_long         | to_name_long           |
|------|-----|-------|----------|---------|-------|-------|------|----------|-------|------------------------|------------------------|
| 14   | 23  | 1.337 | 3.809    | 0.137   | 9.786 | <0.001 | 552  | 0.687    | <0.001 | Prostate Cancer        | Glaucoma               |
| 23   | 14  | 1.056 | 2.875    | 0.138   | 7.631 | <0.001 | 1568 | 0.942    | <0.001 | Glaucoma               | Prostate Cancer        |
| 23   | 29  | 1.691 | 5.425    | 0.485   | 3.488 | <0.001 | 1568 | 0.976    | 0.002 | Glaucoma               | Obesity                |
| 23   | 43  | 0.615 | 1.850    | 0.183   | 3.358 | 0.001 | 1568 | 0.677    | 0.002 | Glaucoma               | Acute Myeloid Leukemia |
| 23   | 44  | 0.901 | 2.463    | 0.265   | 3.404 | 0.001 | 1568 | 0.506    | 0.002 | Glaucoma               | Vitiligo               |
| 23   | 509 | 1.651 | 5.214    | 0.173   | 9.561 | <0.001 | 1568 | 0.279    | <0.001 | Glaucoma               | Med:Progestins         |
| 23   | 63  | 0.769 | 2.158    | 0.239   | 3.212 | 0.001 | 1568 | 0.187    | 0.004 | Glaucoma               | Psoriasis              |
| 29   | 23  | 1.069 | 2.914    | 0.201   | 5.319 | <0.001 | 244  | 0.212    | <0.001 | Obesity                | Glaucoma               |
| 43   | 23  | 0.682 | 1.977    | 0.158   | 4.323 | <0.001 | 350  | 0.069    | <0.001 | Acute Myeloid Leukemia | Glaucoma               |
| 44   | 509 | 0.550 | 1.732    | 0.248   | 2.213 | 0.027 | 424  | 0.473    | 0.057 | Vitiligo               | Med:Progestins         |
| 509  | 23  | 0.529 | 1.697    | 0.102   | 5.179 | <0.001 | 1088 | 0.302    | <0.001 | Med:Progestins         | Glaucoma               |

* Example survival curve: 

<a href="https://github.com/yingeddi2008/scratch/blob/master/14_23.cox.survival.png" target="_blank"><img src="https://github.com/yingeddi2008/scratch/blob/master/14_23.cox.survival.png" alt="example of CoxPH survival curve" style="float: right" width="450"/></a>

In this test example, it shows the probability of getting Glaucoma across time between the population with and without Prostate Cancer. The strata `=0` indicates the population without Prostate Cancer, while strata `=1` indicates with Prostate Cancer. It suggests that people with Prostate Cancer develop Glaucoma faster than those who do not have Prostate Cancer.

* Example residual curve:

<a href="https://github.com/yingeddi2008/scratch/blob/master/14_23.cox.residual.png" target="_blank"><img src="https://github.com/yingeddi2008/scratch/blob/master/14_23.cox.residual.png" alt="example of CoxPH residual curve" style="float: right" width="450"/></a>

For the same Prostate Cancer to Glaucoma disease pair, the residuals lie closely to the diagonal line, meaning the data has a relative good fit for Cox-PH model.

* Example patient drop-off histogram:

<a href="https://github.com/yingeddi2008/scratch/blob/master/14_23.dropoff.png" target="_blank"><img src="https://github.com/yingeddi2008/scratch/blob/master/14_23.dropoff.png" alt="example of dropoff hist" style="float: right" width="450"/></a>

It examines the drop-off subjects distribution for the two populations (having Prostate Cancer and not), and test whether they are similar. Given the P-value = 0.027, which is < 0.05, suggesting that those two distributions are statistically not similar to each other. Users should use caution when considering it to be a meaningful disease correlation.

#### Random Forest survival analysis

* To perform Random Forest survival analysis, run code:

```bash
Rscript a1_analyzer.r --method RF -m /path/to/disease_code_test.tsv --inputfolder /path/to/all_csv_files/ --outputfolder /path/to/all_csv_files/outputs
```

* Example **RF.edge.csv** file output (This won't match the result after you run the code):

| from | to  | Pr_wil | Pr_t  | NonExposedMean | exposedMean | Pr_wil.BH | Pr_t.BH | from_name_long  | to_name_long           |
|------|-----|--------|-------|----------------|-------------|-----------|---------|-----------------|------------------------|
| 14   | 23  | <0.001  | <0.001 | 1501.946       | 982.018     | <0.001     | <0.001   | Prostate Cancer | Glaucoma               |
| 14   | 509 | 0.968  | 0.810 | 1482.072       | 1493.435    | 1.000     | 0.864   | Prostate Cancer | Med:Progestins         |
| 23   | 14  | <0.001  | <0.001 | 1665.489       | 1466.806    | <0.001     | <0.001   | Glaucoma        | Prostate Cancer        |
| 23   | 29  | <0.001  | <0.001 | 1712.080       | 1681.140    | <0.001     | <0.001   | Glaucoma        | Obesity                |
| 23   | 43  | <0.001  | <0.001 | 1661.658       | 1587.812    | <0.001     | <0.001   | Glaucoma        | Acute Myeloid Leukemia |
| 23   | 44  | <0.001  | 0.001 | 1724.847       | 1680.445    | <0.001     | 0.002   | Glaucoma        | Vitiligo               |
| 23   | 509 | <0.001  | <0.001 | 1759.704       | 1557.293    | <0.001     | <0.001   | Glaucoma        | Med:Progestins         |
| 23   | 59  | 0.001  | 0.581 | 1751.787       | 1752.578    | 0.002     | 0.661   | Glaucoma        | Chalazion              |
| 23   | 62  | 0.990  | 0.893 | 1694.242       | 1702.241    | 1.000     | 0.920   | Glaucoma        | Cirrhosis            |
| 23   | 63  | <0.001  | <0.001 | 1653.121       | 1607.370    | <0.001     | <0.001   | Glaucoma        | Psoriasis              |
| 23   | 67  | 0.172  | 0.467 | 1778.002       | 1771.577    | 0.254     | 0.615   | Glaucoma        | Sarcoidosis            |

* Example survival curve: 

<a href="https://github.com/yingeddi2008/scratch/blob/master/14_23.rf.survival.png" target="_blank"><img src="https://github.com/yingeddi2008/scratch/blob/master/14_23.rf.survival.png" alt="example of RF survival curve" style="float: right" width="450"/></a>

* More options:

```
$ Rscript a1_analyzer.r 

 randomForestSRC 2.5.0 
 
 Type rfsrc.news() to see new features, changes, and bug fixes. 
 

       survival randomForestSRC      data.table      doParallel         foreach 
           TRUE            TRUE            TRUE            TRUE            TRUE 
         OIsurv        optparse         ggplot2       ggfortify       gridExtra 
           TRUE            TRUE            TRUE            TRUE            TRUE 
Usage: Rscript a1_analyzer.r -i <inputFileName> -m <metaFileName> [options]

		 >> This is the step 2 of EMR package. << 
	This R script will perform either Cox-PH regression or random forest survival analysis on the valid files in result folder from step 1. 
		Outputs:
		 > CoxPH: an adjacency matrix with all the significant hazard ratios, hazard ratio test p-values, and survival curve and residual graphs.
		 > RF: Wilcoxon and T test p-values, mean time to event in exposed and non-exposed group, and survival curve graphs. 
		*Multiple test correction will be applied to the p-values in both methods.

Options:
	--inputfolder=CHARACTER
		Input folder name. 
		Should be a folder that contains all disease trajectories. Each file is named as NameA_NameB.csv. Default is current directory

	-m CHARACTER, --metafile=CHARACTER
		Meta file for each disease name [required]. 
		First column is the disease ID used in the input file column Disease, second column is the disease name. All the other columns should be additional attribute of the disease, tab separated. Header is needed

	--outputfolder=CHARACTER
		Survival analysis / random forest PNG output folder, default is the <inputfolder>

	--method=CHARACTER
		Survival Analysis method, choises are CoxPH, RF. The default is CoxPH

	-o CHARACTER, --outfile=CHARACTER
		Adjacency matrix output file name, default is <method>.edges.csv

	--duration=INTEGER
		study duration (years), default is 5 years

	-s DOUBLE, --sigcut=DOUBLE
		Significant P-value cutoff. Only p-values that are less than this cutoff will be considered as significant P-values, and kept in the adjacency matrix. Default is 0.05

	-e DOUBLE, --exp_coef=DOUBLE
		Exponentiated coefficients cutoff used in plotting survival curve. Default is 1.0

	-a CHARACTER, --adjust=CHARACTER
		Method of adjusting p-values, choices are holm, hochberg, hommel, bonferroni, BY, and fdr. The default is BH

	-n INTEGER, --ntree=INTEGER
		Number of trees to run random forest, default is 1000

	--pair
		A flag to control for matched subject in CoxPH. Default is false

	-p INTEGER, --processors=INTEGER
		Number of cores/CPUs to use, default is 8

	-h, --help
		Show this help message and exit
```

### Step 3: a2_parse.py
* Parse the CoxPH adjacency matrix into a json object for webpage network display.

```bash
$ python a2_parse.py /path/to/CoxPH.edge.csv
```
* A screenshot of final network display. For interactive network display, go to <a href="http://cbi.lumc.edu/disease/" target="_blank">http://cbi.lumc.edu/disease/</a>.

<a href="https://github.com/yingeddi2008/scratch/blob/master/EMR_display.png" target="_blank"><img src="https://github.com/yingeddi2008/scratch/blob/master/EMR_network.PNG" alt="Screenshot of network disply"/></a>

### Step 4: Copy results to web server directory

* Files: `CoxPH.edge.csv`, `RF.edge.csv`, `all.edges.csv.js`, and all `PNG` files should be in the same folder as the `index.html`.

```bash
cp /path/to/all_csv_files/outputs/* web_server
```

## Version
* Version 1.6.3 Publication version
* Version 1.6 Major improvement on the speed 

## Authors

* Huaiying Lin. M.S., algorithm development, program coding and testing
* Dr. Xiang Gao, theoretical conception and algorithm development
* Dr. Qunfeng Dong, algorithm development
* Kashi Revanna, cytoscape visualization
* Petar Bajic, manuscript drafting and revision
* Michael Zhao, software testing
* Ruichen Rong, algorithm development, program coding and testing

## Error report

Please report any errors or bugs to hlin2@luc.edu.

## License
GNU
