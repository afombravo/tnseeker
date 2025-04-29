<div align="center">

[![Active Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)
![PyPI](https://img.shields.io/pypi/v/tnseeker.svg)
![Docker Image Version (latest by date)](https://img.shields.io/docker/v/afombravo/tnseeker/latest)
![Docker Pulls](https://img.shields.io/docker/pulls/afombravo/tnseeker)

</div>

# Tnseeker
Tnseeker is an advanced pipeline tailored for transposon insertion sequencing (Tn-Seq) analysis. 

It performs an array of tasks: 
1. Read trimming based on the presence of transposon sequences & extraction of associated linked barcodes (if present)
2. Alignment to reference genome (bowtie2)
3. Links genomic locations (using .gff or .gb files as input) with transposon insertion locations
4. Infer essential genes based on global and local (contig wise) transposon insertion densities

What truly distinguishes Tnseeker from other tools is its unique capability to automatically infer and adjust threshold/cutoff parameters. This negates the need for intricate user input, allowing for a more hands-off determination of gene essentiality based on the data. 

Tnseeker is also compatible with any transposon disruption experiment, be it Tn5, HIMAR, or anything else. Tnseeker is thus versatile enough to handle all Tn-Seq datasets.

Tnseeker is under active developement and is available as is. Contact me if you are interested in adapting the program or have any questions. Bugs can be expected. Please report any weird or unintented behaviour. 

## Instalation
There are two main ways of installing tnseeker:


### 1. Recommended installation (using singularity)
1.  Install singularity in your system 

```
conda create -n singularity -c conda-forge singularity -y
```

2.  Active the conda environment
```
conda activate singularity
```

3.  Download the docker image from dockerhub. This will write a singularity container named `tnseeker_latest.sif` into your current work directory.
```bash
singularity pull docker://afombravo/tnseeker:latest
```

4. Start an interactive session of the container. Importantly, you need to `--bind` all the input files. It is easiest if you put all input files into a single folder, as this is how `tnseeker` expects it's input. The results will be written into the same folder.
Note: the `:rw` at the end of the path is crucial for singularity to obtain read/write permission and hence be able to compute.
```
singularity shell --bind /path/to/folder/containing/all/input/files:/input_files:rw \
                  tnseeker_latest.sif
```

5. Start a `tnseeker` run with, for example, the following command:

```
cd /input_files; 
tnseeker \
  --cpu 4 \
  -s ORGANISM_FASTA/GB/GFF_NAME \
  -sd .  \
  -ad .  \
  -at gff \
  -st SE \
  --tn TN_SEQUENCE (ex: AGATTA) \ 
  --m 6 \
  --b \
  --b1 UPSTREAM_BARCODE_SEQUENCE (ex: AGAGA) \
  --b2 DOWNSTREAM_BARCODE_SEQUENCE (ex: ATATAT) \
  --ph 10 \
  --mq 20 \
  --b1m 3 \
  --b2m 3 \
  --ig 100 
```
When using HPC systems it is advisable to include the `--cpu` flag and specify the amount of threads.


---
### 2. Recommended installation #2 (using docker)
1.  Install docker in your system
2.  Download the docker image from dockerhub
```bash
docker pull afombravo/tnseeker:latest
```

3.  Rename to just tnseeker
```bash
docker tag afombravo/tnseeker:latest tnseeker
```

Alternatively, download the docker file from this repo and build it yourself.
```bash
docker build --no-cache -t tnseeker .
```

4.  Start tnseeker docker image with the comand:
```bash
docker run -it -v "<local_path/to/all/your/data>:/data" tnseeker
```

5.  Start tnseeker with:
```bash
tnseeker -sd ./ -ad ./ <ALL OTHER TNSEEKER COMANDS HERE>
```

NOTE: all files required by tnseeker, such as .fasta, .fastq, .gb, or .gff, need to be in the local folder indicated in 4. You then can use the -sd and -ad flags as indicated here in 5.


---
### 3. Alternative installation
The tnseeker pipeline requires Python3, Bowtie2, and BLAST, to be callable from the terminal (and added to path). 

#### For local BLAST
```bash
apt update
apt install ncbi-blast+
```

#### For bowtie2
```bash
apt install bowtie2
```

#### PyPI module
tnseeker can be installed as PyPI module with the folowing:
```bash
pip install tnseeker
```

tnseeker is executable from the command line by typing:

```bash
tnseeker
```

---
## Running Tnseeker

Tnseeker also has a test mode where the blast and Bowtie2 instalations are tested, and a small run on a test dataset is performed.

```bash
tnseeker --tst
```

A quick to use example case is the folowing. See below the meaning of the input arguments:

```bash
tnseeker -s BW25113 -sd ./ -ad ./ -at gb -st SE --tn AGATGTGTATAAGAGACAG --ph 10 --mq 40
```


### File requirements

tnseeker requires several input files:

 1. A '.fastq.gz' file (needs to be .gz)
 
 2. An annotation file in genbank format (.gb), or a .gff (there is an example gff format file in this repo)
 
 3. A FASTA file with the genome under analysis.

---


## Arguments and their meaning:

Flag   | Description
----- | ---
`-h` | Show all parameters and their description.
`-s` | Strain name. Must match the annotation/sequence (FASTA/GFF/GB) file names. TnSeeker will generate a folder with this name where all output will be stored.
`-sd` | The full path to the sequencing files FOLDER.
`--sd_2` | The full path to the pair-ended sequencing files FOLDER (needs to be different from the first folder. In this case keep both ends in different folders).
`-ad` | The full path to the directory with the .gb/.gbk/.gff and .fasta files. TnSeeker searches the directory for file with the name indicated in `-s`
`-at` | Annotation Type (Genbank: `-at gb` / gff: `-at gff`)
`-st` | Sequencing type (Paired-ended (`-st PE`) / Single-ended(`-st SE`)).
`--tst` | Runs TnSeeker in test mode to confirm everything works as expected.
`--tn` | Transposon border sequence (e.g: `--tn GATGTGTATAAGAGACAG`). Required for triming and better/faster mapping. If this parameter is ommited no trimming will be performed and the sequences will be aligned as they are.
`--m` | The number of allowed mismatches in transposon border sequence (default = 0). 
`--k` | If passed, keeps all temporary files (default is disabled: deletes all temporary files).
`--e` | Run only TnSeeker´s essential genes prediction workflow. Requires the "all_insertions_STRAIN.csv" file to have been generated first, and available in the output folder (you can create a folder with the name given in `-s` and place this file there).
`--ne` | Run TnSeeker without the essential genes prediction workflow (outputs the file "all_insertions_STRAIN.csv" file with all annotated transposon insertions, alongside several other quality metrics and plots).
`--t` | Trims to the indicated nucleotides length AFTER finding the transposon sequence. For example: 100 would mean to keep the first 100bp after the transposon (this trimmed read will be used for alignement after).
`--b` | Run with read barcode extraction. Useful when a transposon insertion is linked to a random (or not) barcode in the same read. Enabling this option will output the linkage of all barcode sequences with their respective genome positions and gene annotations.
`--b1` | Upstream barcode delimiting search sequence (example: `--b1 ATC`). This will return any sequence downstream of the input sequence. Needs to be ran together with `--b2` and `--b`.
`--b2` | Downstream barcode delimiting search sequence (example: `--b2 CTA`). This will return any sequence upstream of the input sequence. Together with `--b1` and `--b` these parameters will return the associated transposon barcode.
`--b1m` | Upstream barcode delimiting search sequence mismatches (default is 0).
`--b2m` | Downstream barcode delimiting search sequence mismatches (default is 0).
`--b1p` | Upstream barcode delimiting sequence minimal Phred-score (default is 1).
`--b2p` | Downstream barcode delimiting sequence minimal Phred-score (default is 1).
`--rt` | Read threshold number in absolute reads number (default is 0). All insertions with <= reads will be removed.
`--ph` | Minimal Phred-score (removes reads if any nucleotide has a lower phred-score) (default is Phred-score 1).
`--mq` | Bowtie2 minimal MAPQ threshold (removes reads with a lower MAPQ score) (default is 1).
`--ig` | The number of bp up- and down-stream of any gene that will be ignored before a non maped DNA strech is considered an intergenic region.
`--pv` | Starting p-value threshold for essentiality determination. p-value will be lowered iterativelly based on the optimization of the gold set essential genes. If no gold set genes are detected the value will be the default one (or other if indicated) (Default is 0.05).
`--dut` | The correct essentiality calling of "grey-area" features that have both 'inconclusive' and 'non-essential' sub-gene divisions need special case handling. When the latter two are present, if the number of 'inconclusive domains' > ('inconclusive domains' + 'non-essential') * dut-value (The default is 0.75), then a feature will be demeed 'inconclusive' in its entirity. Otherwise, 'non-essential'.
`--sl5` | 5' gene trimming fraction to be ignored for essentiality determination. For example: `--sl5 0.1` will ignore the insertions in the FIRST 10% of any feature (number between 0 and 1).
`--sl3` | 3' gene trimming fraction to be ignored for essentiality determination. For example: `--sl3 0.9` will ignore the insertions in the LAST 10% of any feature (number between 0 and 1).
`--cpu` | Define the number of threads (must be an integer). Advisable when using HPC systems.

Note: One - denotes mandatory arguments. Two -- indicates options.


---
## Python Dependencies

Tnseeker requires several dependencies, all automatically instalable.
A notable exception is the poibin module, which is available in the current tnseeker folder (you as the user don't need to do anything else), and can be originally be found here: https://github.com/tsakim/poibin


---
### Working modes

tnseeker is composed of 2 submodules: 

1. the initial sequencing processing: Handles the read trimming and alignment, creating a compiled .csv with all found transposon insertions. When individual transposon read associated barcodes are present, these are also extracted. Running only this subpipeline can be achieved by including the `--ne` flag

2. The Essential_finder: Infers gene essentiality from the insertion information found in the previous .csv file. tnseeker can thus be run on a standalone mode if the appropriate .csv and annotation files are indicated. Running only this subpipeline can be achieved by including the `--e` flag
