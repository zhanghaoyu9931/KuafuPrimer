 ![Build](https://img.shields.io/badge/Build-passing-brightgreen) ![Pytorch](https://img.shields.io/badge/Pytorch-V0.0.1-brightgreen) ![license](https://img.shields.io/badge/license-GPL--v3-blue)

## KuafuPrimer: [**Machine learning facilitates the design of 16S rRNA gene primers with minimal bias in bacterial communities**]()

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

We provide two ways to use the KuafuPrimer tool: use repository from GitHub or use the GUI.

###### How to use KuafuPrimer from GitHub

1. Download the program:

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

###### How to use GUI version of KuafuPrimer

TODO YXW: ######

## Quickstart:

###### Design ecosystem-specific primer pair.

User needs to download the [silva.nr_v138_1](https://mothur.org/wiki/silva_reference_files/) dataset for taxonomic classification using mothur tool. These files should be saved as:

1. `Model_data/Silva_ref_data/silva.nr_v138_1.align`
2. `Model_data/Silva_ref_data/silva.nr_v138_1.tax`

Other SILVA dataset files we have processed before need to be downloaded from "Zenodo dataset", and put them in the directory `Model_data/Silva_ref_data/`.The script to produce these files is in the `DataProcess_and_InsilicoEvaluation.ipynb`.

There are some tool requirements before running the program：

1. [mothur](https://mothur.org/) (version == 1.48.0)
2. [muscle](https://drive5.com/muscle) (version == 5.1)
3. [PrimerMatch](https://edwardslab.bmcb.georgetown.edu/software/PrimerMatch.html) (version == 1.0.0)
4. [blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (version == 2.16.0+)
5. [MFEprimer](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (version == 3.2.6)

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

Notably, the input relevant genera profiling file could be generated through the [pipeline](#a-pipeline-for-preprocessing-the-metagenomic-data-of-the-studied-environment.) we used, or provided by users themselves. The designed primers targeting every candidate V-regions will be in the `$demo_output` directory.

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
16sDeepSeg_train_test.ipynb # code for training and evaluation of 16sDeepSeg module
DataProcess_and_InsilicoEvaluation.ipynb # code for preprocessing SILVA dataset and in-silico evaluation
 
```

## Optional function

###### A pipeline for preprocessing the metagenomic data of the studied environment.

Here we privide a pipeline to process metagenomic raw data, please feel free to use your own familiar metagenomic processing workflow instead. This pipeline requires some tools to be installed before:

1. [prinseq-lite](https://github.com/uwb-linux/prinseq) (version == 0.20.4)
2. [KneadData](https://github.com/biobakery/kneaddata#kneaddata-user-manual) (version == 0.12.0)
3. [Usearch](http://www.drive5.com/usearch/) (version == 11.0)
4. [SortMeRNA](https://bioinfo.univ-lille.fr/sortmerna/sortmerna.php) (version == 4.3.6)

`./input/raw_seqs/example_1.fastq` and `./input/raw_seqs/example_2.fastq` are two mates of an example paired-end metagenomic data. And the metagenomic data ids to process should be recorded in `./input/metagenomic_id_list.txt`. To run the pipeline in paired-end mode, run:

```powershell
bash Metagenomic_preprocessing/pipeline_metagenomic.sh \
  ./input \ # directory of metagenomic files
  .R1.raw.fastq.gz .R2.raw.fastq.gz \ # suffix of paired-end data
  ./input/metagenomic_id_list.txt \ # id list

python process_after_pipeline.py \
  --metagenomic_dir ./input \ # directory of metagenomic files
  --id_list ./input/metagenomic_id_list.txt \ # id list
```

This pipeline will output the processed files of the metagenomic data in `./input/` directory. The detailed result files for each sample will be saved in `./input/clean_reads/`. And the integrated abundance table of all samples will be saved in `./input/MetaAbun/`, which will be used as input profiles for the next steps.

###### Use DeepAnno16 to demarcate 16s rRNA gene sequences.

```powershell
# 1. Download the trained 16sDeepSeg modol from [http:***].
# 2. Put the model.pth file in "Model_data/16sDeepSeg_publicated_model/model.pth".
# 3. Run 16sDeepSeg to demarcate the example file:
python 16sDeepSeg/module_output.py -c 16sDeepSeg/config_test.json -r Model_data/16sDeepSeg_publicated_model/model.pth --input input/demo_input_16sDeepSeg.csv --output output/demo_output_16sDeepSeg.csv

```

There will be two files in the output directory: `16sDeepSeg_running_time.txt` records the running time of this procedure, `demo_output_16sDeepSeg.csv` is the 16sDeepSeg output file.

> If you want to design ecosystem specific primer pairs using 16s rRNA gene sequences directly sequenced from the studied ecosystem (which may be absent from the SILVA dataset), you could run the 16sDeepSeg module to demarcate your sequences and afterwards run the design and in-silico PCR procedures

## Citation

TODO: ######

## Contact

If you have any questions, please don't hesitate to ask me: zhanghaoyu9931@pku.edu.cn or hqzhu@pku.edu.cn

## License

The source code of IPEV is distributed as open source under the [GNU GPL v3](https://www.gnu.org/licenses/gpl-3.0.en.html) , the KuafuPrimer program and all datasets  also can be freely available at  [zhulab homepage](https://cqb.pku.edu.cn/zhulab/info/1006/1156.htm)
