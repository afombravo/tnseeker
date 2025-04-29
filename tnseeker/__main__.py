import os
import subprocess
from tnseeker import Essential_Finder,insertions2plots,reads_trimmer,sam2insertions
from tnseeker.extras.helper_functions import cpu,colourful_errors,subprocess_cmd,file_finder
import argparse
from colorama import Fore
import importlib.resources as resources

''' Tnseeker is a pipeline for transposon insertion sequencing (Tn-Seq) analysis. 
    It performs various operations such as read trimming, read
    alignement to a reference genome, random barcode extraction,
    essential gene extraction, and providing several library quality 
    metrics and plots.
    ''' 

def path_finder_seq():
    
    variables["fasta"]=file_finder(variables["annotation_folder"],['*.fasta'],variables["strain"])
    
    variables["directory"] = os.path.join(os.getcwd(), variables["strain"])
    
    if not os.path.isdir(variables["directory"]):
        os.mkdir(variables["directory"])
     
    if variables["sequencing_files"] == None:
        colourful_errors("FATAL",
                 "check that .fastq files exist in the indicated folder.")
        raise Exception
        
    if variables["fasta"] == []:
        colourful_errors("FATAL",
                f"check that the {variables['strain']}.fasta file exist in the indicated folder.")
        raise Exception

    if variables["annotation_type"] == "gb":
        extention = ['*.gb']
        if file_finder(variables["annotation_folder"],extention,variables["strain"]) == []:
            extention = ['*.gbk']
        variables["annotation_file"]=file_finder(variables["annotation_folder"],extention,variables["strain"])
        variables["genome_file"]=file_finder(variables["annotation_folder"],extention,variables["strain"])

    elif variables["annotation_type"] == "gff":
        variables["annotation_file"]=file_finder(variables["annotation_folder"],['*.gff'],variables["strain"])
        variables["genome_file"]=file_finder(variables["annotation_folder"],['*.fasta'],variables["strain"])

    if variables["annotation_file"] == []:
        colourful_errors("FATAL",
            "check that the annotation file exists in the indicated folder.")
        exit()

def bowtie_index_maker():
    
    variables["index_dir"] = f"{variables['directory']}/indexes/"

    if not os.path.isdir(variables["index_dir"]):
        os.mkdir(variables["index_dir"])
            
        send = ["bowtie2-build",
                variables['fasta'],
                f"{variables['index_dir']}{variables['strain']}"]

        subprocess_cmd(send)
    
def bowtie_aligner_maker():
    
    if not os.path.isfile(f'{variables["directory"]}/alignment.sam'):
        if variables["seq_type"] == "SE":
            send = ["bowtie2",
                    "--end-to-end",
                    "-x",f"{variables['index_dir']}{variables['strain']}",
                    "-U",f"{variables['fastq_trimed']}",
                    "-S",f"{variables['directory']}/alignment.sam",
                    "--no-unal",
                    f"--threads {variables['cpus']}",
                    f"2>'{variables['directory']}/bowtie_align_log.log'"]
        else:
            send = ["bowtie2",
                    "--end-to-end",
                    "-x",f"{variables['index_dir']}{variables['strain']}",
                    "-1",f"{variables['fastq_trimed'][0]}",
                    "-2",f"{variables['fastq_trimed'][1]}",
                    "-S",f"{variables['directory']}/alignment.sam",
                    "--no-unal",
                    f"--threads {variables['cpus']}",
                    f"2>'{variables['directory']}/bowtie_align_log.log'"]
        
        subprocess_cmd(send)
        
    else:
        colourful_errors("INFO",
            f"Found {variables['directory']}/alignment.sam, skipping alignment.")

    if variables["remove"]:
        if variables["seq_type"] == "SE":
            os.remove(variables['fastq_trimed'])
        else:
            os.remove(variables['fastq_trimed'][0])
            os.remove(variables['fastq_trimed'][1])
    
def tn_compiler():
    import gzip
    variables["fastq_trimed"] = f'{variables["directory"]}/processed_reads_1.fastq'

    if not os.path.isfile(variables["fastq_trimed"]):
        if file_finder(variables["annotation_folder"],['*.gz'],variables["strain"]) == None:
            
            pathing = file_finder(variables['sequencing_files'],['*.fastq'])
                        
            for file in pathing:
                with open(file, "rb") as firstfile, open(variables["fastq_trimed"],"wb") as secondfile:
                    for line in firstfile:
                        secondfile.write(line)
            
        else:
            pathing = file_finder(variables['sequencing_files'],['*.gz'])
                        
            for file in pathing:
                with gzip.open(file, "rb") as firstfile, open(variables["fastq_trimed"],"wb") as secondfile:
                    for line in firstfile:
                        secondfile.write(line)

def tn_trimmer():
    try:
        variables['sequencing_files']
    except IndexError:
        colourful_errors("FATAL",
            "Make sure that you have selected the correct sequencing type, or that the .gz files are named correctly.")
        exit()

    found_trim = True
    if variables["seq_type"] == "PE":
        variables["fastq_trimed"] = [f'{variables["directory"]}/processed_reads_1.fastq']+\
                                [f'{variables["directory"]}/processed_reads_2.fastq']
        if not (os.path.isfile(variables["fastq_trimed"][0])) & (os.path.isfile(variables["fastq_trimed"][1])):
            found_trim = False
    else:
        variables["fastq_trimed"] = f'{variables["directory"]}/processed_reads_1.fastq'
        if not os.path.isfile(variables["fastq_trimed"]):
            found_trim = False

    if not found_trim:
        reads_trimmer.main(variables["all_variables_path"])
    else:
        colourful_errors("INFO",
            f"Found {variables['fastq_trimed']}, skipping trimming.")

def sam_parser():
    
    if not os.path.isfile(f'{variables["directory"]}/all_insertions_{variables["strain"]}.csv'):
        sam2insertions.main(variables["all_variables_path"])
        
    else:
        colourful_errors("INFO",
            f"Found all_insertions_{variables['strain']}.csv, skipping tn insertion parsing.")

def test_functionalities():
    result_bowtie = subprocess.run(['bowtie2', '-h'], capture_output=True, text=True)
    if result_bowtie.returncode == 0:
        colourful_errors("INFO",
            "Bowtie2 is working as intended.")
    else:
        colourful_errors("FATAL",
            "Bowtie2 is not working as intended. Check instalation and/or that it is on path.")

    result_blast = subprocess.run(['tblastn', '-h'], capture_output=True, text=True)
    if result_blast.returncode == 0:
        colourful_errors("INFO",
            "Blast is working as intended.")
    else:
        colourful_errors("FATAL",
            "Blast is not working as intended. Check instalation and/or that it is on path.")

    if (result_blast.returncode == 0) & (result_bowtie.returncode == 0):
        colourful_errors("INFO",
            "Testing Tnseeker. Please hold, this might take several minutes.")

        data_dir = resources.files('tnseeker').joinpath('data/test/')
        result_full = subprocess.run(["python","-m", "tnseeker", 
                                        "-s","test",
                                        "-sd", data_dir,
                                        "-ad", data_dir,
                                        "-at", "gb",
                                        "-st", "SE",
                                        "--tn", "AGATGTGTATAAGAGACAG",
                                        "--ph", "10",
                                        "--mq", "40",
                                        "--sl5", "0.05", "--sl3", "0.9",
                                        "--k"], capture_output=True, text=True)
        
        if result_full.returncode == 0:
            colourful_errors("INFO",
                "Tnseeker is working as intended.")
        else:
            print(result_full.stdout)
            print(result_full.stderr)
            colourful_errors("FATAL",
                "Tnseeker is not working as intended. Check errors.")

    if (result_blast.returncode == 0) & (result_bowtie.returncode == 0) & (result_full.returncode == 0):
        colourful_errors("INFO",
                " All tests passed.")
    print("\n")
    
def input_parser():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-s",help="Strain name. Must match the annotation (FASTA/GB) file names")
    parser.add_argument("-sd",help="The full path to the sequencing files FOLDER")
    parser.add_argument("--sd_2",help="The full path to the pair ended sequencing files FOLDER (needs to be different from the first folder)")
    parser.add_argument("-ad",help="The full path to the directory with the annotation (FASTA/GB) files")
    parser.add_argument("-at",help="Annotation Type (Genbank)")
    parser.add_argument("-st",help="Sequencing type (Paired-ended (PE)/Single-ended(SE)")
    parser.add_argument("--tn",nargs='?',const=None,help="Transposon border sequence (tn5: GATGTGTATAAGAGACAG). Required for triming and proper mapping")
    parser.add_argument("--m",nargs='?',const=None,help="Mismatches in the transposon border sequence (default is 0)")
    parser.add_argument("--k",nargs='?',const=False,help="Remove intermediate files. Default is yes, remove.")
    parser.add_argument("--e",nargs='?',const=False,help="Run only the essential determing script. required the all_insertions_STRAIN.csv file to have been generated first.")
    parser.add_argument("--t",nargs='?',const=-1,default=-1,help="Trims to the indicated nucleotides length AFTER finding the transposon sequence. For example, 100 would mean to keep the 100bp after the transposon (this trimmed read will be used for alignement after)")
    parser.add_argument("--b",nargs='?',const=False,help="Run with barcode extraction")
    parser.add_argument("--b1",nargs='?',const=False,help="upstream barcode sequence (example: ATC)")
    parser.add_argument("--b2",nargs='?',const=False,help="downstream barcode sequence (example: CTA)")
    parser.add_argument("--b1m",nargs='?',const=False,help="upstream barcode sequence mismatches")
    parser.add_argument("--b2m",nargs='?',const=False,help="downstream barcode sequence mismatches")
    parser.add_argument("--b1p",nargs='?',const=False,help="upstream barcode sequence Phred-score filtering. Default is no filtering")
    parser.add_argument("--b2p",nargs='?',const=False,help="downstream barcode sequence Phred-score filtering. Default is no filtering")
    parser.add_argument("--rt",nargs='?',const=False,help="Read threshold number")
    parser.add_argument("--ne",nargs='?',const=False,help="Run without essential Finding")
    parser.add_argument("--ph", nargs="?", const=1, default=1, help="Phred Score (removes reads where nucleotides have lower phred scores)")
    parser.add_argument("--mq",nargs='?',const=1,default=1, help="Bowtie2 MAPQ threshold")
    parser.add_argument("--ig",nargs='?',const=0,default=0, help="The number of bp up and down stream of any gene to be considered an intergenic region")
    parser.add_argument("--dut",nargs='?',const=0.75,default=0.75, help="Fraction of the minimal amount of 'inconclusive domains' in a gene before the entire gene is deemed too uncertain for essentiality inference")
    parser.add_argument("--pv",nargs='?',const=0.05,default=0.05,help="Essential Finder pvalue threshold for essentiality determination")
    parser.add_argument("--sl5",nargs='?',const=0,default=0,help="5' gene trimming percent for essentiality determination (number between 0 and 1)")
    parser.add_argument("--sl3",nargs='?',const=1,default=1,help="3' gene trimming percent for essentiality determination (number between 0 and 1)")
    parser.add_argument("--tst",nargs='?',const=True,help="Test the program functionalities and dependencies")
    parser.add_argument("--cpu",nargs='?',const=None,help="Define the number of threads (must be and integer)")

    args = parser.parse_args()
                                                             
    print("\n")
    print(f"{Fore.RED} Welcome to{Fore.RESET}")
    print(f"{Fore.RED} ████████╗███╗   ██╗ {Fore.RESET}███████╗███████╗███████╗██╗  ██╗███████╗██████╗ ")
    print(f"{Fore.RED} ╚══██╔══╝████╗  ██║ {Fore.RESET}██╔════╝██╔════╝██╔════╝██║ ██╔╝██╔════╝██╔══██╗")
    print(f"{Fore.RED}    ██║   ██╔██╗ ██║ {Fore.RESET}███████╗█████╗  █████╗  █████╔╝ █████╗  ██████╔╝")
    print(f"{Fore.RED}    ██║   ██║╚██╗██║ {Fore.RESET}╚════██║██╔══╝  ██╔══╝  ██╔═██╗ ██╔══╝  ██╔══██╗")
    print(f"{Fore.RED}    ██║   ██║ ╚████║ {Fore.RESET}███████║███████╗███████╗██║  ██╗███████╗██║  ██║")
    print(f"{Fore.RED}    ╚═╝   ╚═╝  ╚═══╝ {Fore.RESET}╚══════╝╚══════╝╚══════╝╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝")   
    
    variables["version"]="1.2"
    
    print(f"{Fore.RED}            Version: {Fore.RESET}{variables['version']}")
    print("\n")  
    
    if args.tst is not None:
        test_functionalities()
        exit()
        
    if (args.s is None) or (args.sd is None) or (args.ad is None) or (args.at is None) or (args.st is None):
        print(parser.print_usage())
        colourful_errors("FATAL",
                 "No arguments given.")
        raise ValueError

    variables["full"]=True
    if args.e is not None:
        variables["full"] = False

    variables["trim"]=False
    variables["tn_mismatches"] = 0 
    if args.tn is not None:
        variables["trim"] = True

        if args.m is not None:
            variables["tn_mismatches"] = int(args.m)   

    variables["remove"]=True
    if args.k is False:
        variables["remove"]=False

    variables["trimmed_after_tn"]=args.t
        
    variables["barcode"]=False
    if args.b is not None:
        variables["barcode"] = True

    variables["intergenic_size_cutoff"]=0
    if args.ig is not None:
        variables["intergenic_size_cutoff"] = int(args.ig)
    
    variables["barcode_up"] = None
    variables["barcode_down"] = None
    variables["barcode_up_miss"] = 0 
    variables["barcode_down_miss"] = 0 
    variables["barcode_up_phred"] = 1
    variables["barcode_down_phred"] = 1
    if variables["barcode"]:

        if args.b1p is not None:
            variables["barcode_up_phred"] = int(args.b2p)
            if variables["barcode_up_phred"] < 1:
                variables["barcode_up_phred"]  = 1
                
        if args.b2p is not None:
            variables["barcode_down_phred"] = int(args.b2p)
            if variables["barcode_down_phred"] < 1:
                variables["barcode_down_phred"]  = 1
        
        if args.b1m is not None:
            variables["barcode_up_miss"] = int(args.b1m)

        if args.b2m is not None:
            variables["barcode_down_miss"] = int(args.b2m)
        
        if args.b1 is not None:
            variables["barcode_up"] = str(args.b1).upper()

        if args.b2 is not None:
            variables["barcode_down"] = str(args.b2).upper()

    variables["read_threshold"]=False
    variables["read_value"] = 0
    if args.rt is not None:
        variables["read_threshold"] = True
        variables["read_value"] = args.rt
        
    variables["essential_find"]=True
    if args.ne is not None:
        variables["essential_find"]=False

    variables["pvalue"]=float(args.pv)
    variables["domain_uncertain_threshold"]=float(args.dut)
    variables["subdomain_length_up"]=float(args.sl5)
    variables["subdomain_length_down"]=float(args.sl3)

    if args.cpu is not None:
        variables["cpus"]=int(args.cpu)
    else:
        variables["cpus"]=cpu()
    
    if args.st == "PE":
        variables["sequencing_files_r"] = args.sd_2
        
    variables["sequencing_files"] = args.sd
    variables['phred']=int(args.ph)
    variables['MAPQ']=int(args.mq)
    variables["strain"]=args.s
    variables["annotation_type"]=args.at.lower()
    variables["annotation_folder"]=args.ad
    variables["seq_type"]=args.st
    variables["sequence"]=args.tn
    variables["cmd_used"] = " ".join(f"--{key}" if isinstance(value, bool) and value else f"--{key} {value}" for key, value in vars(args).items() if value is not None)

def variables_initializer():

    global variables
    variables = {}
    input_parser()
    path_finder_seq()
    variables["all_variables_path"] = os.path.join(variables['directory'],"cmd_input.txt")
    
    print(f"{Fore.YELLOW} -- Parameters -- {Fore.RESET}\n")
    with open(variables["all_variables_path"],'w+') as current:
        for key in variables:
            current.write(str(key)+' : '+str(variables[key])+'\n')
            if ("all_variables_path" not in key) & ("cmd_used" not in key):
                print(f"{Fore.GREEN} {key}:{Fore.RESET} {variables[key]}")
    print(f"\n{Fore.YELLOW} ---- {Fore.RESET}\n")

def final_report_compiler():
    final_report_path = os.path.join(variables['directory'],"final_report.txt")

    trim_report = file_finder(variables["directory"],["*.log"],"trimming_log")
    if trim_report != []:
        with open(trim_report, "rb") as firstfile, open(final_report_path,"ab+") as secondfile:
            for line in firstfile:
                secondfile.write(line)
        os.remove(trim_report)

    align_report = file_finder(variables["directory"],["*.log"],"bowtie_align_log")
    if align_report != []:
        copy_permit = False
        with open(align_report, "rb") as firstfile, open(final_report_path,"ab+") as secondfile:
            for line in firstfile:
                if "reads; of these:" in str(line,"utf-8"):
                    copy_permit = True
                if copy_permit:
                    secondfile.write(line)
        os.remove(align_report)

    quality_report = file_finder(variables["directory"],["*.log"],"library_stats_")
    if quality_report != []:
        with open(quality_report, "rb") as firstfile, open(final_report_path,"ab+") as secondfile:
            for line in firstfile:
                secondfile.write(line)
        os.remove(quality_report)

    essentiality_report = file_finder(variables["directory"],["*.log"],"Essentiality_stats")
    if essentiality_report != []:
        with open(essentiality_report, "rb") as firstfile, open(final_report_path,"ab+") as secondfile:
            for line in firstfile:
                secondfile.write(line)
        os.remove(essentiality_report)

def main():

    variables_initializer()
    
    if variables["full"]:
        bowtie_index_maker()

        if variables["trim"]:
            colourful_errors("INFO",
                "Getting that .fastq ready.")
            tn_trimmer()

        else:
            colourful_errors("INFO",
                "Compiling those .fastq.")
            if variables["seq_type"] == "SE": #directly skipping transposon based triming currently only works for SE
                tn_compiler()
        
        colourful_errors("INFO","Aligning reads to the reference genome.")
        bowtie_aligner_maker()

        sam_parser()
        if variables["remove"]:
            os.remove(f"{variables['directory']}/alignment.sam")
        insertions2plots.main(variables["all_variables_path"])
    
    if variables["essential_find"]:
        colourful_errors("INFO","Infering essential genes.")
        Essential_Finder.main(variables["all_variables_path"])

    final_report_compiler()
    
    colourful_errors("INFO","Analysis Finished.")

if __name__ == "__main__":
    main()