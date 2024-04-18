# Tnseeker
Tnseeker is an advanced pipeline tailored for transposon insertion sequencing (Tn-Seq) analysis. It performs an array of tasks: from read trimming and alignment to associating genomic locations with transposon insertions and inferring essential genes based on transposon insertion densities. Additionally, Tnseeker is adept at extracting barcodes from raw fastq files and linking them to corresponding transposon genomic locations for subsequent analysis. What truly distinguishes Tnseeker from other tools is its unique capability to automatically infer and adjust threshold/cutoff parameters. This negates the need for intricate user input, allowing for a more precise determination of gene essentiality based on the data. Compatible with any transposon disruption experiment, Tnseeker efficiently mitigates transposon-specific biases, including those seen with HIMAR. Hence, Tnseeker is versatile enough to handle all Tn-Seq datasets.

Tnseeker is under active developement and is available as is. Contact me if you are interested in using the program or have any questions. Bugs can be expected. Please report any weird or unintented behaviour. 

## Requirements
The tnseeker pipeline requires Python3, Bowtie2, and BLAST, to be callable from the terminal (and added to path). 
Tnseeker also requires several python dependencies, all automatically installed. A notable exception is the poibin module, which is available in the current tnseeker folder (you as the user don't need to do anything), and can be originally be found here: https://github.com/tsakim/poibin

## Instalation
`sudo apt update`

#### For local BLAST
`sudo apt install ncbi-blast+`

#### For bowtie2
`sudo apt install bowtie2=2.4.4-1`

#### Tnseeker 
Tnseeker requires Python>=3.9 installed.
Tnseeker can then be installed as a PyPI module with the folowing:

`pip install tnseeker`

## Executing 

tnseeker is executable from the command line by typing:

`python -m tnseeker`

An example use case is the folowing. See below the meaning of the input arguments:

`python -m tnseeker -s BW25113 -sd '/your/data/directory/folder_with_fastq.gz_files' -ad /your/annotations/directory/ -at gb -st SE --tn AGATGTGTATAAGAGACAG --ph 10 --mq 40`

## Optional Arguments:

  -h, --help   show this help message and exit

  -s S         Strain name. Must match the annotation (FASTA/GB) file
               names

  -sd SD       The full path to the sequencing files FOLDER

  --sd_2 SD_2  The full path to the pair ended sequencing files FOLDER (needs
               to be different from the first folder)

  -ad AD       The full path to the directory with the .gb and .fasta files

  -at AT       Annotation Type (Genbank)

  -st ST       Sequencing type (Paired-ended (PE)/Single-ended(SE)

  --tn [TN]    Transposon border sequence (tn5: GATGTGTATAAGAGACAG). Required for triming and proper mapping

  --m [M]      Mismatches in the transposon border sequence (default is 0)

  --k [K]      Remove intermediate files. Default is yes, remove.

  --e [E]      Run only the essential determing script. required the
               all_insertions_STRAIN.csv file to have been generated first.

  --t [T]      Trims to the indicated nucleotides length AFTER finding the
               transposon sequence. For example, 100 would mean to keep the
               100bp after the transposon (this trimmed read will be used for
               alignement after)

  --b [B]      Run with barcode extraction

  --b1 [B1]    upstream barcode sequence (example: ATC)

  --b2 [B2]    downstream barcode sequence (example: CTA)

  --b1m [B1M]  upstream barcode sequence mismatches

  --b2m [B2M]  downstream barcode sequence mismatches

  --b1p [B1P]  upstream barcode sequence Phred-score filtering. Default is no
               filtering

  --b2p [B2P]  downstream barcode sequence Phred-score filtering. Default is
               no filtering
  --rt [RT]    Read threshold number

  --ne [NE]    Run without essential Finding

  --ph [PH]    Phred Score (removes reads where nucleotides have lower phred
               scores)

  --mq [MQ]    Bowtie2 MAPQ threshold

  --ig [IG]    The number of bp up and down stream of any gene to be
               considered an intergenic region

  --pv [PV]    Essential Finder pvalue threshold for essentiality
               determination

  --dut [DUT]  fraction of the minimal amount of 'too small domains' in a gene before the entire gene is deemed
               uncertain for essentiality inference
  
  --sl5 [SL5]  5' gene trimming percent for essentiality determination (number
               between 0 and 1)

  --sl3 [SL3]  3' gene trimming percent for essentiality determination (number
               between 0 and 1)
               

### File requirements

tnseeker requires several input files:

 1. A '.fastq.gz' file (needs to be .gz)
 
 2. An annotation file in genbank format (.gb)
 
 3. A FASTA file with the genome under analysis (needs to be .fasta).


### Working modes

tnseeker is composed of 2 submodules: 

1. the initial sequencing processing: Handles the read trimming and alignment, creating a compiled .csv with all found transposon insertions. When individual transposon read associated barcodes are present, these are also extracted.

2. The Essential_finder: Infers gene essentiality from the insertion information found in the previous .csv file. tnseeker can thus be run on a standalone mode if the appropriate .csv and annotation files are indicated. 
