# High-resolution immune atlas of murine healthy prostate and prostate cancer from neoplasia to adenocarcinoma stage in the inducible PTEN(i)pe-/- mice

## Article information

**Title:** High-resolution immune atlas of murine healthy prostate and prostate cancer from neoplasia to adenocarcinoma stage in the inducible PTEN(i)pe-/- mice

**Authors:** 
Despoina Pervizou1, Joanna De Chiara1, Lionel Spinelli1, Lionel Chasson1, Frédéric Fiore, Marc Bajénoff1, Bernard Malissen1, Daniel Metzger, Gilles Laverny and Sandrine Henri1&

1 Aix-Marseille Univ, Centre National de la Recherche Scientifique (CNRS), Institut National de la Santé et de la Recherche Médicale (INSERM), Centre d'Immunologie de Marseille-Luminy (CIML), Marseille, France.

& Corresponding author (e-mail: sandrine.henri@inserm.fr)

**Summary:**
To establish a high-resolution immune atlas of healthy prostate and prostate cancer we performed scRNA-seq on sorted tissue-resident CD45+ immune cells from the inducible adult PTEN(i)pe-/- mice at the Prostatic Intraepithelial Neoplasia stage (PIN) and Adenocarcinoma stages, and age-matched healthy PTENL2/L2 control mice at 3 months and 9 months post-tamoxifen administration respectively. A few minutes prior sacrifice, mice were injected intravenously with fluorescent anti-CD45.2 antibodies to perform blood/tissue partitioning, allowing to exclude cells present in the vasculature and ensure to study solely of the immune cells present in prostate parenchyma. Cells were extracted from DLV prostate lobes by enzymatic digestion. Prostate samples from three PTEN(i)pe−/− mice at 3 months post-tamoxifen were pooled prior to sorting. Similarly, prostate samples from five  sex-matched PTENL2/L2 control littermates at 3 months post-tamoxifen were also pooled prior to sorting. The same approach was applied for the 9 months samples of PTEN(i)pe−/− and control mice respectively. Tissue resident CD45+ immune cells were sorted by fluorescence-activated cell sorting (FACS). As we found previously that neutrophils were prevalent within the leukocyte infiltrate from PIN to adenocarcinoma stages of Pten(i)pe−/−prostate, neutrophils were sorted separately from the other CD45+ immune cells and cells were mixed again at a ratio of 25% neutrophils to 75% CD45+ non-neutrophil immune cells prior 10X Genomics Chromium System encapsulation and scRNA sequencing. We included 8 independent samples encompassing 2 replicates of PtenL2/L2 control at 3 months, PtenL2/L2 control at 9 months, Pten(i)pe−/− at 3 months PIN and Pten(i)pe−/− at 9 months adenocarcinoma. 

**DOI:**

---

---

## Goal of the GitHub
This GitHub project contains the instructions and materials to reproduce the analyses reported in the article (and more).
Source code (scripts and dockerfile) are available in the GitHub repository. Processed data, analysis results, and built Docker/Singularity images are available for download from Recherche Data Gouv. Instructions to reproduce the analyses are provided below.

---

## Description of the datasets

As described in the article, there are 8 libraries in this study, each of which has been sequenced twice, generating 16 datasets.

The libraries are the result of the composition of two genotypes (PTENL2/L2 and PTEN(i)pe−/−), two timepoints ( 3 months and 9 months Tamoxifen) and two replicates (1 and 2).

A dataset name contains the information required to identify it uniquely as follows:
```
[FLOWCELL_NAME]_[CONDITION][TIMEPOINT]_[REPLICATE]
```

where:

 * FLOWCELL_NAME is one of "230127_VH00228_144_AAAGNCJHV" and "230601_VH00228_181_AAAGNF7HV"
 * CONDITION is one of "CTRL" (PTENL2/L2) and "PTEN" (PTEN(i)pe−/−)
 * TIMEPOINT is one of "3mo" (3 mo Tamoxifen) and "9mo" (9 mo Tamoxifen)
 * REPLICATE is one of "1" and "2"

| Dataset name | title | genotype	| treatment	| description
| :--------------- |:--------------- |:-------|:-------|:-------|
| 230127_VH00228_144_AAAGNCJHV_CTRL3mo_1 | Mus Musculus Prostate Control 3 months post tam Replicate 1 Sequencing Run 1 | PTENL2/L2 | 3 mo Tamoxifen | WT prostate (mouse age: 5 months)
| 230127_VH00228_144_AAAGNCJHV_CTRL3mo_2	| Mus Musculus Prostate Control 3 months post tam  Replicate 2 Sequencing Run 1 | PTENL2/L2 | 3 mo Tamoxifen | WT prostate (mouse age: 5 months)
| 230127_VH00228_144_AAAGNCJHV_CTRL9mo_1	| Mus Musculus Prostate Control 9 months post tam Replicate 1 Sequencing Run 1 | PTENL2/L2 | 9 mo Tamoxifen | WT prostate (mouse age: 11 months)
| 230127_VH00228_144_AAAGNCJHV_CTRL9mo_2	| Mus Musculus Prostate  Control 9 months post tam Replicate 2 Sequencing Run 1 | PTENL2/L2 | 9 mo Tamoxifen | WT prostate (mouse age: 11 months)
| 230127_VH00228_144_AAAGNCJHV_PTEN3mo_1	| Mus Musculus Prostate PIN 3 months post tam Replicate 1 Sequencing Run 1 | PTEN(i)pe−/− | 3 mo Tamoxifen | Tumor prostate at PIN stage (mouse age: 5 months)
| 230127_VH00228_144_AAAGNCJHV_PTEN3mo_2	| Mus Musculus Prostate PIN 3 months post tam Replicate 2 Sequencing Run 1 | PTEN(i)pe−/− | 3 mo Tamoxifen | Tumor prostate at PIN stage (mouse age: 5 months)
| 230127_VH00228_144_AAAGNCJHV_PTEN9mo_1	| Mus Musculus Prostate Adenocarcinoma 9 months post tam Replicate 1 Sequencing Run 1 | PTEN(i)pe−/− | 9 mo Tamoxifen | Tumor prostate at Adenocarcinoma stage (mouse age: 11 months)
| 230127_VH00228_144_AAAGNCJHV_PTEN9mo_2	| Mus Musculus Prostate Adenocarcinoma 9 months post tam Replicate 2 Sequencing Run 1 | PTEN(i)pe−/− | 9 mo Tamoxifen | Tumor prostate at Adenocarcinoma stage (mouse age: 11 months)
| 230601_VH00228_181_AAAGNF7HV_CTRL3mo_1	| Mus Musculus Prostate Control 3 months post tam Replicate 1 Sequencing Run 2 | PTENL2/L2 | 3 mo Tamoxifen | WT prostate (mouse age: 5 months)
| 230601_VH00228_181_AAAGNF7HV_CTRL3mo_2	| Mus Musculus Prostate Control 3 months post tam  Replicate 2 Sequencing Run 2 | PTENL2/L2 | 3 mo Tamoxifen | WT prostate (mouse age: 5 months)
| 230601_VH00228_181_AAAGNF7HV_CTRL9mo_1	| Mus Musculus Prostate Control 9 months post tam Replicate 1 Sequencing Run 2 | PTENL2/L2 | 9 mo Tamoxifen | WT prostate (mouse age: 11 months)
| 230601_VH00228_181_AAAGNF7HV_CTRL9mo_2	| Mus Musculus Prostate  Control 9 months post tam Replicate 2 Sequencing Run 2 | PTENL2/L2 | 9 mo Tamoxifen | WT prostate (mouse age: 11 months)
| 230601_VH00228_181_AAAGNF7HV_PTEN3mo_1	| Mus Musculus Prostate PIN 3 months post tam Replicate 1 Sequencing Run 2 | PTEN(i)pe−/− | 3 mo Tamoxifen | Tumor prostate at PIN stage (mouse age: 5 months)
| 230601_VH00228_181_AAAGNF7HV_PTEN3mo_2	| Mus Musculus Prostate PIN 3 months post tam Replicate 2 Sequencing Run 2 | PTEN(i)pe−/− | 3 mo Tamoxifen | Tumor prostate at PIN stage (mouse age: 5 months)
| 230601_VH00228_181_AAAGNF7HV_PTEN9mo_1	| Mus Musculus Prostate Adenocarcinoma 9 months post tam Replicate 1 Sequencing Run 2 | PTEN(i)pe−/− | 9 mo Tamoxifen | Tumor prostate at Adenocarcinoma stage (mouse age: 11 months)
| 230601_VH00228_181_AAAGNF7HV_PTEN9mo_2	| Mus Musculus Prostate Adenocarcinoma 9 months post tam Replicate 2 Sequencing Run 2 | PTEN(i)pe−/− | 9 mo Tamoxifen | Tumor prostate at Adenocarcinoma stage (mouse age: 11 months)


---

## Description of the main analysis

During the analysis process of the datasets, several steps were applied to combine and merge the datasets. Here is a quick description.

### First step : Individual dataset analysis

All the datasets were analyzed separately to validate their quality and get a first insight into the cell heterogeneity. 

The analysis code and analysis results for this step are in folders with the dataset name.

### Second step : Aggregate resequenced datasets

 For each dataset, all the FASTQ files from both sequencing runs were used to run a new Cell Ranger count analysis. These aggregated data were analyzed separately to take advantage of the increased read saturation.

The analysis code and analysis results for this step are in folders with a name carrying the two flowcell names. 

For instance,  __VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL3mo_1__ is the folder name for the aggregation of the two datasets __230127_VH00228_144_AAAGNCJHV_CTRL3mo_1__ and __230601_VH00228_181_AAAGNF7HV_CTRL3mo_1__.

### Third step : Merge dataset by timepoints

The aggregated datasets from same time points were merged together to compare conditions at each time points. Thanks to the tightly controlled experimental setup, a simple data merge was enough to align the data (no integration process was required).

The analysis code and analysis results for this step are in folders with a name carrying the two flowcells names and the two conditions names with "Merge" as prefix. For instance __Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo__ if the merge of __VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL3mo_1__, __VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL3mo_2__ and __VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN3mo_1__, __VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN3mo_2__

### Fourth step : Merge all datasets

All the aggregated datasets were merged together to compare both time points and condition effects.

The analysis code and analysis results for this step are in a folder called __Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo__.

---

## Structure of the data and code

### Folders of the datasets

The project is organized by dataset. Each dataset has its own folder. Due to the analysis steps described above, when downloading the code and data (see below), you will obtain several sub-folders with names as below:

```
.
└── CD45pos_Prostate_OCT2022
    ├── 230127_VH00228_144_AAAGNCJHV_CTRL3mo_1
    ├── 230127_VH00228_144_AAAGNCJHV_CTRL3mo_2
    ├── 230127_VH00228_144_AAAGNCJHV_CTRL9mo_1
    ├── 230127_VH00228_144_AAAGNCJHV_CTRL9mo_2
    ├── 230127_VH00228_144_AAAGNCJHV_PTEN3mo_1
    ├── 230127_VH00228_144_AAAGNCJHV_PTEN3mo_2
    ├── 230127_VH00228_144_AAAGNCJHV_PTEN9mo_1
    ├── 230127_VH00228_144_AAAGNCJHV_PTEN9mo_2
    ├── 230601_VH00228_181_AAAGNF7HV_CTRL3mo_1
    ├── 230601_VH00228_181_AAAGNF7HV_CTRL3mo_2
    ├── 230601_VH00228_181_AAAGNF7HV_CTRL9mo_1
    ├── 230601_VH00228_181_AAAGNF7HV_CTRL9mo_2
    ├── 230601_VH00228_181_AAAGNF7HV_PTEN3mo_1
    ├── 230601_VH00228_181_AAAGNF7HV_PTEN3mo_2
    ├── 230601_VH00228_181_AAAGNF7HV_PTEN9mo_1
    ├── 230601_VH00228_181_AAAGNF7HV_PTEN9mo_2
    ├── Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo
    ├── Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo
    ├── Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_9mo
    ├── VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL3mo_1
    ├── VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL3mo_2
    ├── VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL9mo_1
    ├── VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL9mo_2
    ├── VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN3mo_1
    ├── VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN3mo_2
    ├── VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN9mo_1
    └── VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN9mo_2
``` 

### Folders inside the datasets

Each dataset folder contains an ordered series of sub-folders. Theses sub-folders are organized in two sets : the code and the data and analysis results.

 * From the GitHub, you will be able to download the code set of files. It contains the analysis code (__03_Script__) and the Docker container definition file (__02_Container__). Some dataset also contains a Snakemake workflow definition (__04_Worflow__). For instance : 

```
.
└── CD45pos_Prostate_OCT2022
    ├── Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo
    │   ├── 02_Container
    │   ├── 03_Script
    │   ├── 04_Workflow
```

 * From Recherche Data Gouv, you will be able to download the data and analysis results set of files. It contains the reference files (_01_Reference__), the compiled docker images (_02_Container__) and the analysis results (05_Output)
 
```
.
└── CD45pos_Prostate_OCT2022
    ├── Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo
    │   ├── 01_Reference
    │   ├── 02_Container
    │   ├── 05_Output
```

---
---

## Prepare the environments

In order to prepare the environment for analysis execution, it is required to:

- Clone the GitHub repository and set the WORKING_DIR environment variable
- Download the data
- Install the Docker environment to run the analysis interactively
- (optional) Install Snakemake and Singularity to run the analysis workflow

Below you will find detailed instruction for each of these steps.

---

#### Clone the GitHub repository

Use you favorite method to clone this repository in a chosen folder. This will create a folder **CD45pos_Prostate_OCT2022** with all the source code.

#### Set the Working Dir variable

Then, you must set an environment variable called **WORKING_DIR** with a value set to the path to this folder.

For instance, if you have chosen to clone the Git repository in __"/home/spinellil/workspace"__, then the **WORKING_DIR** variable will be set to __"/home/spinellil/workspace/CD45pos_Prostate_OCT2022"__

**On linux:**

```
    export WORKING_DIR=/home/spinellil/workspace/CD45pos_Prostate_OCT2022
```

#### Add your working dir in the code

The code uses variables that are stored in different "parameters" file. One important variable is the PATH_PROJECT which indicate to the code where your project is stored.
You have to modify this variable in the code to reflect your project setup. Each dataset has a file called **globalParams.R** in the subfolder **03_Script**

For instance:

```
    CD45pos_Prostate_OCT2022
    │
    ├── 230127_VH00228_144_AAAGNCJHV_CTRL3mo_1
    │   │
    │   └── 03_Script/globalParams.R

```

Edit those files and in each of them locate the line defining the **PATH_PROJECT** variable and change its value to the same value as the **WORKING_DIR** variable you defined before. Then save the files.

```
PATH_PROJECT = "/home/spinellil/workspace/CD45pos_Prostate_OCT2022"
```

---

### Download the data

The raw FASTQ files are available on GEO. The rest of the data are available on Recherche Data Gouv into 4 datasets :

- 230127_VH00228_144_AAAGNCJHV_AllSamples (https://doi.org/10.57745/NVZQ25) : All samples from first sequencing
- 230601_VH00228_181_AAAGNF7HV_AllSamples (https://doi.org/10.57745/NVZQ25) : All samples from second sequencing
- VH00228_144_AAAGNCJHV_181_AAAGNF7HV_AllAggregations (https://doi.org/10.57745/R1VNP3) : All samples from aggregation of the two sequencing
- VH00228_144_AAAGNCJHV_181_AAAGNF7HV_AllMerges (https://doi.org/10.57745/4EFY2Z) : Merges of samples by Conditions and by Conditions and Timepoints

You can download all of them or part of then according to your needs. Download the folders into WORKING_DIR, in order to keep the correct folder structure.

#### Download the analysis result

The analysis results can be found in the folder called "05_Ouput" in each dataset folder. These subfolders contains a series of subfolders, one for each analysis step. The same subfolders can be found in the 03_Script subfolders, since the analysis code and the analysis output have a bijective relation. 

**Note :** the analysis results contain also the Cell Ranger pre-processing results (including the count matrix and the bam files).

For instance:

```
    CD45pos_Prostate_OCT2022
    │
    ├── Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo$
        └── 03_Script
            ├── 01_DatasetIntegration
            ├── 02_CompareConditions
            ├── 03a_AnalyzeCellTypes_Bcell
            ├── 03b_AnalyzeCellTypes_Tcell
            ├── 03c_AnalyzeCellTypes_Neutrophil
            ├── 03d_AnalyzeCellTypes_Myeloid
            ├── 04a_CompareConditions_Bcell
            ├── 04b_CompareConditions_Tcell
            ├── 04c_CompareConditions_Neutrophil
            ├── 04d_CompareConditions_Myeloid
            ├── 05d_GlobalHeterogeneity_SciGeneX_Myeloid
            ├── 06d_CompareConditions_SciGeneX_Myeloid
            ├── 07d_DynamicAnalysis_Velocyto_Myeloid
```

#### Download the container images

The container images (tar.gz file for Docker and sif file for Singularity) can be found in the 02_Container sub-folder of each dataset folder. Singulairty images can be use directly while docker images must be loaded on your system (see below)

---

### Install the Docker environment to run the analysis interactively

#### Install Docker

You need install Docker on your system to take advantage of interactive analysis environment with Rstudio, follow the instructions here : https://docs.docker.com/get-docker/


#### Load the images

Once done, locate the tar.gz files of the docker images and use the following command

```
docker load -i <image_name>.tar.gz
```

#### Run the analysis individually using Docker

If you have loaded the docker images (see above), you can use Rstudio in Docker to run the analysis individually.

To start a docker container, use the following command:

```
docker run -d -p 8787:8787 -v /$WORKING_DIR:/$WORKING_DIR -e PASSWORD=<PASSWORD> -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g)  <IMAGE_NAME>
```

where:

* <PASSWORD> is a simple string you will use as password to login into Rstudio
* <IMAGE_NAME> is the Docker image name to run

One started, you can open an internet browser and use the URL https://127.0.0.1:8787.

At the login prompt, enter the name of the user session you are connected with and the password you type in place of <PASSWORD>. You are now in a Rstudio environment and the container is able to connect to the **WORKING_DIR**
of your system. Inside you will find the project files. To tun the analysis, look at the scripts "launch_report_compilation.R" that is the main entry point to run the analysis.

### (optional) Install Snakemake and Singularity to run the analysis workflow

To run the analysis using the Snakemake workflow (available only for some datasets), you need to:

 * install Snakemake to run the complete analysis workflow. Use your preferred solution : https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

 * install Singularity v3 on your system to run the complete analysis using the Snakemake Workflow. Follow the instructions here : https://sylabs.io/guides/2.6/admin-guide/

Then use the snakemake executable and the snakefile in the "04_Workflow" sub-folder to run the complete analysis.





