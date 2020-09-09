# Welcome to SilhouetteRank

This toolkit contains silhouetteRank, a flexible method for finding spatially variable genes. It computes a score based on silhouette coefficient of binarized gene expression data. It allows users to specify multiple spatial widths and integrate them in a Fisher's test. 

silhouetteRank is written in Python 3.

## Installation

### Prerequisites
* eva R package [1]
* EnvStats R package
* qvalue R package
* GNU parallel (UNIX)

[1] eva package needs to be installed through file:
```
install.packages("eva_0.2.5.tar.gz", source="file", repos="NULL")
```
Download the file eva_0.2.5.tar.gz [here](https://cran.r-project.org/src/contrib/Archive/eva/eva_0.2.5.tar.gz)

### Install
```bash
pip3 install --user --no-deps --no-cache --upgrade silhouetteRank
```
This will install the Python package locally to user's directory under `~/.local/lib`.


## Usage

There are two ways of using silhouetteRank.
One way is to run on **single-machine** using multiple CPU cores. For this, the GNU parallel library is used to achieve parallelization.

The other way is to run through **SLURM scheduler**. By going this route, we can distribute jobs across multiple nodes in a computing cluster. This is suitable when there are a large number (~100) of computing nodes available.

## Single-machine

Example command:
```bash
~/.local/bin/silhouette_rank_one -x expression.txt -c Xcen.good -r 0.95 0.99 -e 0.005 0.01 0.05 0.10 0.30 -m dissim -p 8 -a /usr/bin -o results -v -q 10
```

**Table 1**. Explanations of `silhouette_rank_one`: 
| param | explanation |
| ----- | ------------------ |
| -x | Expression file |
| -c | Cell position file |
| -r | Array. Float. Local spatial distance weighting constant. (recommend 0.95 - 0.995) |
| -e | Array. Float. Top proportion of cells to binarize to 1 (0 - 1) |
| -m | dissim or sim. Use dissimilarity matrix (default) or similarity matrix |
| -p | Integer. Number of CPU cores to be used for GNU parallel |
| -a | Path to GNU parallel binary (e.g. /usr/bin) |
| -q | Chunk size (advanced). Leave as default of 10. |
| -o | Output directory |
| -v | Flag. Verbose |

Example of expression file [here](http://spatialdataset.com). Example of cell position file [here](http://spatialgiotto.com)

It is a good idea to look into `$output/logs` directory to track the progress of the program.

The final result is saved in the file `$output/silhouette.overall.pval.txt`.

## SLURM scheduler

The SLURM method is even faster than single-machine as it can utilize lots of computing machines to distribute the computation.

```bash
python3 ~/.local/lib/python3.6/site-packages/silhouetteRank/prep.py -x expression.txt -c Xcen.good -r 0.95 0.99 -e 0.005 0.01 0.05 0.10 0.30 -m dissim -o results -v
```
The explanations of parameters are much the same as the `silhouette_rank_one` shown in above Table 1.
This program `prep.py` prepares the necessary directory structure and cache files in the output directory for main steps of the program.

Go to the output directory (specified by `-o` in `prep.py`).
```bash
cd results
```
Create a SLURM **job submission script** `do_one.sh`. This script computes null distribution of silhouette scores. 
```bash
vim do_one.sh
```
Press `i` (for insert). Paste the following content:

```bash
#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=16000                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=bernardzhu@gmail.com   # Email to which notifications will be sent

module load R/3.6.1
rbp_p=$1
top=$2
q_index=$3
echo $1 $2 $3
~/.local/bin/silhouette_rank_random -r $rbp_p -e $top -m dissim -o . -q $q_index -v
```
Press `ESC`. Then `:wq` to save and quit vim. 

Create a **second** job submission script `do_one_real.sh`. This script calculates silhouette scores for actual dataset.
```bash
vim do_one_real.sh
```
Press `i` (for insert), and paste the following content:
```bash
#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=16000                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=bernardzhu@gmail.com   # Email to which notifications will be sent

module load R/3.6.1
rbp_p=$1
top=$2
echo $1 $2
~/.local/bin/silhouette_rank_main -x ../expression.txt -c ../Xcen.good -e $top -r $rbp_p -m dissim -o . -v
```
Press `ESC`. Then `:wq` to save and quit vim. 

Note that in the job submission script. You can adjust the time, memory usage, and the queue name according to your needs. The time (5hr) and resources (16GB memory) shown work for a dataset with 25,000 cells, 11,000 genes. 

Submit all the jobs:
```bash
for i in 0.95 0.99; do for j in 0.005 0.01 0.05 0.1 0.3; do for k in `seq 0 9`; do sbatch ./do_one.sh $i $j $k; done; done; done
```
```bash
for i in 0.95 0.99; do for j in 0.005 0.01 0.05 0.1 0.3; do sbatch ./do_one_real.sh $i $j; done; done
```
In this example, we request one job per combination of (i, j, k). Hence there will be 100 jobs running. In the case of `do_one_real.sh`, we request one job per (i, j) for a total of 10 jobs.

Once all the jobs are complete, combine the results together to generate a final score list.
Note this step does not need to run using SLURM scheduler.
```bash
python3 ~/.local/lib/python3.6/site-packages/silhouetteRank/combine_2.py -i . -v -r 0.95 0.99 -e 0.005 0.01 0.05 0.1 0.3 -m dissim 
```

**Table 2**. Explanations of `combine.py`: 
| param | explanation |
| ----- | ------------------ |
| -i | The results directory (containing directories like result_5000_0.95_0.300) |
| -v | Verbose |
| -r | Array. Float. Local spatial distance weighting constant. (recommend 0.95 - 0.995) |
| -e | Array. Float. Top proportion of cells to binarize to 1 (0 - 1) |
| -m | dissim or sim. Use dissimilarity matrix (default) or similarity matrix |

The final result is saved in the file `silhouette.overall.pval.txt` located in the results directory.

Note: the parameter settings for `-r` and `-e` work well for most datasets. Users can safely use the same commands for all datasets.

