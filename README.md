<img src="https://user-images.githubusercontent.com/52743495/173834089-526c540a-df4b-452f-964e-26104bf6f261.png" width="350" />

**JRC_seeker** is a tool for identifying jointly regulated CpGs (JRCs) in the human methylome. These regions can be classified into a number of genomic phenomenon, such as imprinting regions or methylation quantitative trait loci (mQTLs). Developed by the Genetic Identification Lab at Erasmus Medical Center, a BAM file of WGBS reads can be inputted and the result of this tool is a list of JRC locations and their associated p-values. JRC_seeker is built using a Snakemake pipeline that combines Python scripts with Linux shell commands.

## Requirements

```
Operating system: tested on Ubuntu 18.04.6 LTS (Bionic Beaver)
R: tested on R version 3.6.1 (2019-07-05) -- "Action of the Toes"
Python: Python 3.9.12
RAM requirements: Do not attempt to run without at least 50 GB of RAM.
Runtime: Approximately 1 day for 98 GB BAM file and 3.1 GB reference genome, using 20 cores.
```

## Dependencies

### [conda](https://www.anaconda.com/products/individual)

Download and install conda if you do not have it already on your machine.

```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
bash Anaconda3-2021.11-Linux-x86_64.sh
```

### [mamba](https://github.com/mamba-org/mamba)

Install Mamba into your Conda-based python distribution

```
conda install -n base -c conda-forge mamba
```

Activate the Conda base environment (which now includes Mamba).

```
conda activate base
```

### [snakemake](https://snakemake.readthedocs.io/) (at least v4.3.1)

Create a new conda environment called `jrc_seeker` with snakemake and python 3.9 in it.

```
mamba create -c conda-forge -c bioconda -n jrc_seeker snakemake python=3.9
```

**Option:**
If conda-forge did not work for you, simply create a conda environment like this and use the pip command options below:

```
conda create -n jrc_seeker python=3.9 snakemake
```

Activate the `jrc_seeker` conda environment.

```
conda activate jrc_seeker
```

Check whether Snakemake is succesfully installed by running the following command:

```
snakemake --help
```

### [biopython-1.79](https://biopython.org/docs/1.79/api/Bio.html)

Install the following packages: biopython=1.79 .

```
conda install -c conda-forge biopython=1.79
```

or

```
pip install "biopython==1.79"
```

### [tabix](https://github.com/samtools/htslib)

```
conda install -c bioconda tabix
```

or

```
pip install tabix
```

### [pandas](https://pandas.pydata.org/)

```
conda install -c anaconda pandas
```

or

```
pip install pandas
```

### [bgzip](https://github.com/xbrianh/bgzip.git)

```
conda install -c bioconda bgzip
```

or

```
pip install bgzip
```

### [ChromHMM](http://compbio.mit.edu/ChromHMM/)

Quick instructions on running ChromHMM:

1. Install Java 1.5 or later if not already installed.
2. Download and unzip the ChromHMM.zip file using the following code snippit:

```
wget http://compbio.mit.edu/ChromHMM/ChromHMM.zip
unzip ChromHMM.zip
```

### [bedtools](https://bedtools.readthedocs.io/en/latest/)

```
conda install -c bioconda bedtools
```

or

```
pip install bedtools
```

### [samtools](http://www.htslib.org/doc/samtools.html)

Recommended version: 1.14

```
conda install -c bioconda samtools
```

or

```
pip install "samtools==1.14"
```

### [samblaster](https://github.com/GregoryFaust/samblaster)

```
conda install -c bioconda samblaster
```

or

```
pip install samblaster
```

### [BISCUIT](https://huishenlab.github.io/biscuit/)

```
conda install -c bioconda biscuit
```

### R

If not installed already, be sure you have R version 4.1.2 or greater.

### R dependencies

To install R dependencies, open R by running:

```bash
R
```

And run the following R commands:

```r
# To install
install.packages('data.table')
install.packages('parallel')
install.packages('MASS')

# To verify that installation was succesful
library(data.table)
library(parallel)
library(MASS)
```

## Downloading JRC_Seeker

Clone the repository with the following command:

```
mkdir jrc_seeker
cd jrc_seeker
git clone https://github.com/b-kolar/jrc_seeker.git
```

## Process Overview

1. **Generate methylation data** - _pileup, meth_info, format_meth_data_
2. **Binarize methylation data** - _binarize_
3. **ChromHMM genome segmentation** - _learnmodel_
4. **BinPolish** - _binpolish_assets, label_states, binpolish_

### Generate methylation data

BISCUIT is used to create a pileup VCF of DNA methylation and genetic information. After generating the VCF file, with both genetic and methylation information, beta values and coverage are extracted using the `vcf2bed` command to study the methylation levels at sequenced CpGs. The output of this step is a BED file of the format below, where the columns represent chromosome, CpG start base, CpG end base, beta methylation value (proportion of reads methylated), and total read coverage.

```
chr1    10469   10470   0.625   8
chr1    10471   10472   0.444   9
chr1    10484   10485   0.889   9
chr1    10489   10490   1.000   10
chr1    10493   10494   0.875   8
```

### Binarize methylation data

From the above BED file, the counts of methylated and unmethylated reads at every CpG location are calculated using the `format_methylation.py` script. Using the `format_methylation.py` script CpGs counts are combined using 200-bp bins and a binarized track is formed per chromosome for all chromosomes that contain methylation data:

```
methylated  unmethylated
0           1
0           0
1           1
1           0
1           1
```

The first column represents if the bin contains methylated reads `(0 = False, 1 = True)` and the second column represents if the bin contains unmethylated reads `(0 = False, 1 = True)`. Thus: `0   1` is an unmethylated bin, `1   0` is a methylated bin, and `1   1` is an intermediately methylated bin. This bins are binarized from count data using their methylation value. A methylation value less than 0.2 is unmethylated and above 0.8 is methylated. Between these values is intermediately methylated. These binarization thresholds can be changed in the `config.json` file using the `lower_im_methylation_bound` and the `upper_im_methylation_bound` variables. The bin size can also be adjusted in the `config.json` file using the `bin_size` variable. ChromHMM, the genome segmentation software used in the following step recommends a bin size of 200. Changing this bin size here will also change it for the ChromHMM step.

A threshold is also set to set bins with low coverage to a no-data state, aka `0   0`. If there are less that 3 methylated/unmethylated counts (not CpGs), then this bin is set to the low-coverage state. This threshold can be set by changing the `k_threshold` variable in the `config.json` file.

### ChromHMM genome segmentation

ChromHMM, a genome segmentation and annotation software tool, is leveraged to identify intermediately methylated regions within the provided dataset. Ultimately, this step summarizes the binarized bins into larger regions throughout the genome.

Finding intermediately methylated regions is necessary for finding JRCs. ChromHMM uses a multivariate hidden Markov Model to infer states from binarized data, allowing these regions to be identified and for the rest of the dataset/genome to be removed (as this is not of interest for the purpose of finding CpGs). ChromHMM segments all chromosomes provided into state regions and annotates them with a number. Each of these numbers corresponds to one of the four states of interest to us (no-data, unmethylated, methylated, intermediately methylated).

![image](https://user-images.githubusercontent.com/52743495/174035708-41f2e402-666b-48dd-845a-892f3d7194d9.png)

The output of ChromHMM is a list of segments that correspond to a given state. Using the emission matrix outputted by ChromHMM, the states can be annotated and the regions corresponding to intermediately methylated regions can be identified.

### BinPolish

After the ChromHMM segmentation, there are often a large number of intermediately methylated regions, some of which are very small (200 bp). The goal of BinPolish is to remove some of these sparse, small regions or to summarize them into larger intermediately methylated blocks, thus reducing the large number of small intermediately methylated regions into a set of larger, robust ones. Essentially, we are polishing intermediately methylated segments.

![BinPolish](https://user-images.githubusercontent.com/52743495/174038503-de253ac7-e8c3-4e08-86a5-7d6f0d60af16.png)

Intermediately methylated regions are ignored that overlap with known regions that have low mappability (source) or are documented blacklisted regions (source). The number of CpGs per region is also calculated using the original data from BISCUIT. Following this, a the main steps employed that merge or remove intermediately methylated (IM) regions are:

1. Merge IM regions that are separated by 200bp
2. If small IM regions are within two larger regions, classify as the state of the earlier large region
3. If two IM regions are separated by a region that has less than one CpG, merge these IM regions
4. If IM region is classified wrongly, re-classify
5. Merge IM regions that are separated by 200bp (again)
6. If a state <=600bp is within two IM states that have higher cpg density, merge
7. States with CpG density <= 2 are turned to no-data states
8. If small IM regions are within two larger regions, classify as the state of the earlier large region (again)
9. Remove small regions (<= 200 bp)
#   S a m p l e - C o d e - - - J R C _ S e e k e r  
 