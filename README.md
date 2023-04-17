# tnseeker
Tnseeker is a pipeline for transposon insertion sequencing (Tn-Seq) analysis. It performs various operations such as read trimming, alignement, and essential gene inference.

## Installation
Tnseeker currently exist as an unpublished PyPI module. You can download the installation wheel from the current repository ('tnseeker-1.0.0.tar.gz') and install it in a Linux environment by typing: 

`pip install "/your/tnseeker/directory/tnseeker-1.0.0.tar.gz" `
 
## Requirements
The tnseeker pipeline requires both Python3 and Bowtie2 to be callable from the terminal. 

## Executing 
tnseeker is executable from the command line by typing:

`python -m tnseeker`

An example use case is the folowing. See below the meaning of the input arguments:

`python -m tnseeker -s BW25113 -sd '/your/data/directory/file.fastq.gz' -ad /your/annotations/directory/ -at gb -st SE --tn AGATGTGTATAAGAGACAG --ph 40 --mq 10 --b`

## Optional Arguments:

  -h, --help   show this help message and exit
  
  -s S         Strain name. Must match the annotation (FASTA/GFF/GB) file
               names
               
  -sd SD       The full path to the sequencing file
  
  --sd_2 SD_2  The full path to the pair ended sequencing file
  
  -ad AD       The full path to the directory with the annotation
               (FASTA/GFF/GB) files
               
  -at AT       Annotation Type (Genbank/gff)
  
  -st ST       Sequencing type (Paired-ended (PE)/Single-ended(SE)
  
  --tn [TN]    Transposon border sequence (Himar: 'ACTTATCAGCCAACCTGT'; tn5:
               'GATGTGTATAAGAGACAG'). Required for triming and proper mapping
               
  --k [K]      Remove intermediate files. Default is yes, remove.
  
  --e [E]      Run only the essential determing script. required the
               all_insertions_STRAIN.csv file to have been generated first.
               
  --t [T]      Run without searching for a transposon sequence, and triming
               the .fastq file
               
  --b [B]      Run with barcode extraction
  
  --rt [RT]    Read threshold number
  
  --ne [NE]    Run without essential Finding
  
  --ph [PH]    Phred Score (removes reads where nucleotides have lower phred
               scores)
               
  --mq [MQ]    Bowtie2 MAPQ threshold
  
  --pv [PV]    Essential Finder pvalue threshold for essentiality
               determination
               
  --sl5 [SL5]  5' gene trimming percent for essentiality determination (number
               between 0 and 1)
               
  --sl3 [SL3]  3' gene trimming percent for essentiality determination (number
               between 0 and 1)
               
## Dependencies

tnseeker requires several dependences, all instalable via `pip` commands, except the notable exception of the poibin module, which is available in the current tnseeker folder, and can be originally be found here: https://github.com/tsakim/poibin
