 ![Build](https://img.shields.io/badge/Build-passing-brightgreen) ![Pytorch](https://img.shields.io/badge/Pytorch-V0.0.1-brightgreen) ![license](https://img.shields.io/badge/license-GPL--v3-blue)

## KuafuPrimer: Machine learning facilitates the design of 16S rRNA gene primers with minimal bias in bacterial communities

![0](./pic/logo.png)

## Contents

- [Introduction](#Introduction)
- [Installation](#Installation)
- [Quickstart](#Quickstart)
- [Citation](#Citation)
- [Contact](#Contact)
- [License](#License)

## Introduction

KuafuPrimer is a machine learning-aided method that learns community characteristics from several samples to design 16S rRNA gene primers with minimal bias for microbial communities. It is built on **Python 3.9.0**,  **Pytorch 1.12.0**.

## Installation

File tree:

```
.
├── 16sDeepSeg/
├── input/
├── Metagenomic_process/
├── Model_data/
├── output/
├── README.md
├── requirements_conda.txt
├── requirements_pip.txt
├── classify_seqs.sh
├── KuafuPrimer.py
├── KuafuPrimerFuncs.py
├── Insilico_eva_primers.py
└── Screen_best_PP.py
```

###### How to install KuafuPrimer from GitHub

1. Clone the repository:

   ```powershell
   git clone https://github.com/zhanghaoyu9931/KuafuPrimer.git
   ```
2. You can configure the operating environment (please create new environment to avoid unnecessary trouble) by using:

   ```powershell
   1. conda create -n my_env_name python=3.9
   2. source activate my_env_name
   3. pip install -r requirements_pip.txt (or conda install --yes --file requirements_conda.txt)
   ```

   Notably, you may need to install [**torch==1.12.0+cu113**](https://pytorch.org/get-started/previous-versions/) manually.

###### Dependencies

User needs to download the [silva.nr_v138_1](https://mothur.org/wiki/silva_reference_files/) dataset for taxonomic classification using mothur tool. These files should be saved as:

1. `Model_data/Silva_ref_data/silva.nr_v138_1.align`
2. `Model_data/Silva_ref_data/silva.nr_v138_1.tax`

Other dataset files the KuafuPrimer requires are provided in the Zenodo (10.5281/zenodo.13829178), and user needs to download and put them in the directory `Model_data/Silva_ref_data/`. The script to produce these files is in the `DataProcess_and_InsilicoEvaluation.ipynb`.

There are some tool requirements before running the KuafuPrimer：

1. [mothur](https://mothur.org/) (version == 1.48.0)
2. [muscle](https://drive5.com/muscle) (version == 5.1)
3. [PrimerMatch](https://edwardslab.bmcb.georgetown.edu/software/PrimerMatch.html) (version == 1.0.0)
4. [blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (version == 2.16.0+)
5. [MFEprimer](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (version == 3.2.6)

## Quickstart:

###### Design ecosystem-specific primer pair.

To design the ecosystem specific primer pairs from the relevant genera profiling of the studied environment, you need to provide a csv file containing the relevant genera profiling as `demo_input_KuafuPrimer_ge_profile.csv` and run:

```
# set the input and output path
demo_input=input/relevant_genera_profile_input/demo_input_EcoPrimer_ge_profile.csv
demo_output=output/demo_output_EcoPrimer_ge_profile

# run the program
python KuafuPrimer.py \
    --input $demo_input \ # Input file
    --out_root $demo_output \ # Output dir
    --input_type 'genera_profiling' \ # Input data type
    --target_vs 'v1v2;v1v3;v3v4;v4v4;v4v5;v6v8;v5v6;v5v7;v7v8' \ # The target V-regions.
    --extend_bp_num 50 \ # Number of extensive bp that will be included to design primer.
    --num_every_spe 5 \ # Number of selected seqs for every genus.
    --NGS_mode Single_end \ # The metagenomic data type (pair-end or single-end).
```

Notably, the input relevant genera profiling file could be generated through the [pipeline](#Optional-functions) we used, or provided by users themselves. The designed primers targeting every candidate V-regions will be in the `$demo_output` directory.

###### In-silico PCR of the designed primer pairs and screen for the primer with minimal bias for the studied communities.

To evaluate the perfromance of designed primer pairs by in-silico PCR and screen out the optimal primer pair for the studied ecosystem, run:

```
# set the parameters
demo_PCR_output=output/demo_PCR_output
K_num=3

# run the in-silico PCR program
python Insilico_eva_primers.py \
    --envi_forEva $demo_output \ # The samples in which the primer performance is to be evaluated.
    --primers_forEva $demo_output \ # The designed primers to be evaluated.
    --target_vs 'v1v2;v1v3;v3v4;v4v4;v4v5;v6v8;v5v6;v5v7;v7v8' \ # The target V-regions.
    --K $K_num \ # The permitted mismatch numbers.
    --output $demo_PCR_output \ # The output dir.

# parse the in-silico PCR res and screen out the best primer pair
python Screen_best_PP.py \
    --pcr_dir $demo_PCR_output'_K'$K_num \
    --rk_by Genus_accuracy

```

The input for `envi_forEva` and `primers_forEva` is the output directory of primer design  procedure. This program will output the designed communities specific primer pair, the selected V-region and in-silico PCR accuracy of it:

```
Designed ecosystem specific primer pair for output/demo_PCR_output is ['GYCACAYTGGRACTGAGA', 'GGACTACCAGGGTATCTAA'], targeting v3v4 V-region, with 97.27% in-silico PCR amplicon accuracy.
```

The detailed in-silico PCR performance and meta-information of every condidate primer pair are recorded in `detail_Genus_accuracy.csv` and `pri_metainfo.csv` files within the `$demo_PCR_output'_K'$K_num` directory.

## Reproduce results in the publication

Here we provide the scripts for reproducing the results of our paper:

```
DeepAnno16_train_test.ipynb # code for training and evaluation of DeepAnno16 module
DataProcess_and_InsilicoEvaluation.ipynb # code for preprocessing SILVA dataset and in-silico evaluation
 
```

## Optional functions

###### A pipeline for preprocessing the metagenomic data of the studied environment.

Here we privide a pipeline to process metagenomic raw data, please feel free to use your own familiar metagenomic processing workflow instead. This pipeline requires some tools to be installed before:

1. [sra-tools](https://github.com/ncbi/sra-tools) (version == 3.0.2)
2. [fastp](https://github.com/OpenGene/fastp) (version == 0.23.4)
3. [bowtie2](https://github.com/BenLangmead/bowtie2) (version == 2.3.4.1)
4. [Kraken2](https://github.com/DerrickWood/kraken2) (version == 2.1.3)
5. [Bracken](https://github.com/jenniferlu717/Bracken) (version == 2.9)

And some external databases need to be downloaded:

1. Kraken databases
   ```powershell
   cd ./Metagenomic_preprocessing/

   ## build standard kraken2 database
   DBNAME=./database/ncbi_db_kraken2
   kraken2-build --standard --db $DBNAME --threads 32
   kraken2-build --download-taxonomy --db $DBNAME --threads 32
   kraken2-build --download-library bacteria --db $DBNAME --threads 32
   kraken2-build --build --db $DBNAME --threads 32

   # bracken
   KRAKEN_DB=$DBNAME
   READ_LEN=150
   bracken-build -d ${KRAKEN_DB} -t 32 -l ${READ_LEN}

   ## build silva kraken2 database
   DBNAME=./database/silva_db_kraken2
   kraken2-build --special silva --db $DBNAME --threads 32
   kraken2-build --download-taxonomy --db $DBNAME --threads 32
   kraken2-build --download-library bacteria --db $DBNAME --threads 32
   kraken2-build --build --db $DBNAME --threads 32

   # bracken
   KRAKEN_DB=$DBNAME
   READ_LEN=150
   bracken-build -d ${KRAKEN_DB} -t 32 -l ${READ_LEN}

   ```
2. Human reference database for bowtie2
   ```powershell
   cd ./Metagenomic_preprocessing/

   wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip # download the zip file
   unzip GRCh38_noalt_as.zip -d ./database/GRCh38/ # unzip the file
   ## users may need to rename the decompressed files corresponding to the names used in the 'pipeline_metagenomic.sh'.
   ```

`./input/raw_seqs/example_1.fastq` and `./input/raw_seqs/example_2.fastq` are two mates of an example paired-end metagenomic data. And the metagenomic data ids to process should be recorded in `./input/metagenomic_id_list.txt`. To run the pipeline in paired-end mode, run:

```powershell
cd ./Metagenomic_preprocessing/

bash pipeline_metagenomic.sh \
  ../input \ # directory of metagenomic files
  .R1.raw.fastq.gz .R2.raw.fastq.gz \ # suffix of paired-end data
  ../input/metagenomic_id_list.txt \ # id list

python process_after_pipeline.py \
  --metagenomic_dir ../input \ # directory of metagenomic files
  --id_list ../input/metagenomic_id_list.txt \ # id list
```

This pipeline will output the processed files of the metagenomic data in `./input/` directory. The detailed result files for each sample will be saved in `./input/clean_reads/`. And the integrated abundance matrix of all samples will be saved in `./input/MetaAbun/`, which will be used as input profiles for the next steps. There will be two types of abundance matrixes named as `merged_abun_table_ncbi.tsv` (classified by kraken with ncbi database) and `merged_abun_table.tsv` (classified by kraken with silva database), **and we recommend to use the `merged_abun_table_ncbi.tsv`.**

###### Use DeepAnno16 to demarcate 16S rRNA gene sequences.

```powershell
# 1. Download the trained DeepAnno16 modol from Zenodo (10.5281/zenodo.13829178).
# 2. Put the model.pth file in "Model_data/DeepAnno16_publicated_model/model.pth".
# 3. Run DeepAnno16 to demarcate the example file:
python DeepAnno16/module_output.py -c DeepAnno16/config_test.json -r Model_data/DeepAnno16_publicated_model/model.pth --input input/demo_input_DeepAnno16.csv --output output/demo_output_DeepAnno16.csv

```

There will be two files in the output directory: `DeepAnno16_running_time.txt` records the running time of this procedure, `demo_output.csv` is the DeepAnno16 output file.

> If you want to design primer pairs using 16s rRNA gene sequences directly sequenced from the studied ecosystem (which may be absent from the SILVA dataset), you could run the DeepAnno16 module to demarcate your sequences and afterwards run the design and in-silico PCR procedures **(TODO)**

## Citation

KuafuPrimer: Machine learning facilitates the design of 16S rRNA gene primers with minimal bias in bacterial communities

## Contact

If you have any questions, please don't hesitate to ask me: zhanghaoyu9931@pku.edu.cn or hqzhu@pku.edu.cn
