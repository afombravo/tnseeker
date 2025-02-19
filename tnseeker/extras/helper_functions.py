import multiprocessing
from datetime import datetime
from colorama import Fore
import csv
import subprocess
import os
import glob
from regex import findall

def cpu():
    """
    Determine the number of available CPU cores, leaving one free.

    Parameters
    ----------
    None

    Returns
    -------
    int
        Number of CPU cores available for processing.
    """
    c = multiprocessing.cpu_count()
    if c >= 2:
        c -= 1
    return c


def colourful_errors(warning_type, error):
    """
    Print a color-coded error message with a timestamp.

    Parameters
    ----------
    warning_type : str
        Type of warning (e.g., "FATAL", "WARNING").
    error : str
        Error message to be displayed.

    Returns
    -------
    None
    """
    warning_colour = Fore.GREEN
    if warning_type == "FATAL":
        warning_colour = Fore.RED
    elif warning_type == "WARNING":
        warning_colour = Fore.YELLOW

    print(f"{Fore.BLUE} {datetime.now().strftime('%c')}{Fore.RESET} [{warning_colour}{warning_type}{Fore.RESET}] {error}")

def csv_writer(output_file_path,output_file):
    """
    Write a list of rows to a CSV file.

    Parameters
    ----------
    output_file_path : str
        Path to the output CSV file.
    output_file : list
        List of rows (each row as a list) to be written.

    Returns
    -------
    None
    """
    with open(output_file_path, "w", newline='') as output:
        writer = csv.writer(output)
        writer.writerows(output_file)

def subprocess_cmd(command):
    """
    Run a command using subprocess and return the output.

    Parameters
    ----------
    command : list
        List of command-line arguments.

    Returns
    -------
    str
        Output of the command if successful, or error message if failed.
    """
    try:
        return subprocess.check_output(command)
    except subprocess.CalledProcessError as e:
        return e.output.decode()
    
def file_finder(folder, ext, search_term = None):
    """
    Find files in a given folder with specific extensions and an optional search term.

    Parameters
    ----------
    folder : str
        Path to the folder to search.
    ext : list
        List of file extensions to look for.
    search_term : str, optional
        Substring to filter filenames (default is None).

    Returns
    -------
    list or str
        List of matching file paths or a single file path if search_term is specified.
    """
    pathing = []
    for exten in ext:
        for filename in glob.glob(os.path.join(folder, exten)):
            if search_term is None:
                pathing.append(filename) 
            else:
                find = filename.find(search_term) 
                if find != -1: 
                    return filename
    return pathing

def variables_parser(var_file):
    """
    Parse a file of "key : value" entries and return a dictionary of key-value pairs.

    Parameters
    ----------
    var_file : str
        Path to the variable file.

    Returns
    -------
    dict
        Dictionary containing parsed key-value pairs.
    """
    variables = {}
    with open(var_file) as current:
        for line in current:
            key,value = line[:-1].split(" : ")
            variables[key] = value
    return variables

def adjust_spines(ax, spines,x,y):
    """
    Adjust the spines of a matplotlib axis for better visibility.

    Parameters
    ----------
    ax : matplotlib axis
        The axis to modify.
    spines : list
        List of spines to adjust (e.g., ["left", "bottom"]).
    x : tuple
        Bounds for the x-axis spine.
    y : tuple
        Bounds for the y-axis spine.

    Returns
    -------
    None
    """
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 2.5))  
            if (loc == 'left') or (loc == 'right'):
                spine.set_bounds(y)
            else:
                spine.set_bounds(x)
        else:
            spine.set_color('none') 

def gb_parser(file_line):
    """
    Parse GenBank file information into a structured dictionary.

    Parameters
    ----------
    file_line : BioPython SeqFeature
        A feature line from a GenBank file.

    Returns
    -------
    dict
        Dictionary with parsed genomic feature details.
    """
    output = {}
    output["start"] = file_line.location.start
    output["end"] = file_line.location.end
    output["domain_size"] = output["end"] - output["start"]
    output["orientation"] = file_line.location.strand
    
    try:
        output["ID"] = file_line.qualifiers['locus_tag'][0]
    except KeyError:
        output["ID"] = None

    if output["orientation"] == 1:
        output["orientation"] = "+"
    else:
        output["orientation"] = "-"
    
    try:
        if 'product' in file_line.qualifiers:
            output["product"] = file_line.qualifiers['product'][0]
        else:
            output["product"] = file_line.qualifiers['note'][0]
    except KeyError:
        output["product"] = None
        
    for key, val in file_line.qualifiers.items():   
        if "pseudogene" in key:
            output["gene"] = output["ID"]
            break #avoids continuing the iteration and passing to another key, which would make "gene" assume another value
        elif "gene" in key:
            output["gene"] = file_line.qualifiers['gene'][0]
            break #avoids continuing the iteration and passing to another key, which would make "gene" assume another value
        else:
            output["gene"] = output["ID"]

    return output

def gff_parser(file_line):
    """
    Parse GFF file information into a structured dictionary.

    Parameters
    ----------
    file_line : list
        A line from a GFF file split into fields.

    Returns
    -------
    dict
        Dictionary with parsed genomic feature details.
    """
    output = {}

    if "##FASTA" in file_line[0]:
        return

    if "#" not in file_line[0][:3]: #ignores headers

        output["start"] = int(file_line[3])
        output["end"] = int(file_line[4])
        output["domain_size"] = output["end"] - output["start"]
        
        features = file_line[8].split(";") #gene annotation file

        feature = {}
        for entry in features:
            entry = entry.split("=")
            if len(entry) == 2:
                feature[entry[0]] = entry[1].replace("\n","")
        if "gene" in feature:
            output["gene"]=feature["gene"]
        if "Name" in feature:
            output["gene"]=feature["Name"]
        else:
            output["gene"]=feature["ID"]

        output["ID"] = feature["ID"]
        
        output['product'] = None
        if 'product' in feature:
            output["product"] = feature["product"]
            
        if output["gene"] == ".":
            output["gene"] = f'{file_line[2]}_{output["start"]}_{output["end"]}'
        
        output["contig"] = file_line[0]
        output["orientation"] = file_line[6] #orientation of the gene
    
    return output

def inter_gene_annotater(genes,contigs,intergenic_size_cutoff):
    """
    Annotate intergenic regions based on gene positions and contig lengths.

    Parameters
    ----------
    genes : list
        List of gene metadata tupples as follows (start,end,orientation,gene,contig)
    contigs : dict
        Dictionary of contig lengths.
    intergenic_size_cutoff : int
        Cutoff for defining intergenic regions.

    Returns
    -------
    dict
        Dictionary of annotated intergenic regions.
    """

    genes = list(dict.fromkeys(genes))
    genes.sort(key=lambda x: (x[-1], x[0])) #sort by start position of the gene and contig

    ir_annotation = {}
    count = 0
    for i,gene in enumerate(genes[:-1]):
        contig = gene[-1]
        if contig == genes[i+1][-1]: #same contigs
            gene_up_start_border = intergenic_size_cutoff + gene[1]
            gene_down_start_border = genes[i+1][0] - intergenic_size_cutoff
            domain_size = gene_down_start_border - gene_up_start_border
            if domain_size >= 1:
                count += 1
                ir_annotation[f'IR_{count}_{gene[3]}_{gene[2]}_UNTIL_{genes[i+1][3]}_{genes[i+1][2]}'] = (gene_up_start_border,gene_down_start_border,domain_size,contig)

        if contig != genes[i+1][-1]:

            circle_closer = gene[1] + intergenic_size_cutoff 
            domain_size = contigs[contig] - circle_closer
            if domain_size >= 1:
                count += 1
                ir_name = f'IR_{count}_{gene[3]}_{gene[2]}_contig_{contig}_-end'
                if ir_name not in ir_annotation:
                    ir_annotation[ir_name] = (circle_closer,contigs[contig],domain_size,contig)

            domain_size = genes[i+1][0] - intergenic_size_cutoff
            if domain_size  >= 1:
                count += 1
                ir_name = f'IR_{count}_contig_{genes[i+1][-1]}_-start_{genes[i+1][3]}_{genes[i+1][2]}'
                if ir_name not in ir_annotation:
                    ir_annotation[ir_name] = (0,domain_size,domain_size,genes[i+1][-1])

    circle_closer = genes[-1][1] + intergenic_size_cutoff 
    domain_size = contigs[contig] - circle_closer
    if domain_size >= 1:
        count += 1
        ir_name = f'IR_{count}_{genes[-1][3]}_{genes[-1][2]}_contig_{contig}_-end'
        if ir_name not in ir_annotation:
            ir_annotation[ir_name] = (circle_closer,contigs[contig],domain_size,contig)

    domain_size = genes[0][0] - intergenic_size_cutoff
    if domain_size  >= 1:
        count += 1
        ir_name = f'IR_{count}_contig_{contig}_-start_{genes[0][3]}_{genes[0][2]}'
        if ir_name not in ir_annotation:
           ir_annotation[ir_name] = (0,domain_size,domain_size,contig)
    
    return ir_annotation

def bowtie2parser(sam,mapq,flag_list):
    """
    Parse a Bowtie2 SAM alignment file and extract relevant fields.

    Parameters
    ----------
    sam : list
        A line from a SAM file split into fields.
    mapq : int
        Minimum mapping quality threshold.
    flag_list : list
        List of flag values to filter reads.

    Returns
    -------
    dict
        Dictionary with parsed alignment information.
    """
    output = {}
    if (sam[0][0] != "@") and (sam[2] != '*'): #ignores headers and unaligned contigs
        output["local"] = sam[3]
        output["sequence"] = sam[9] 
        output["contig"] = sam[2]
        output["flag"] = int(sam[1])
        output["cigar"] = sam[5]
        output["map_quality"] = float(sam[4])
        output["border"], output["orientation"] = "",""
        output["aligned_reads"] = 1
        multi = "XS:i:" in sam #multiple alignemnts

        if (output["flag"] in flag_list) & (multi==False) & (output["map_quality"] >= mapq): #only returns aligned reads witht he proper flag score
            
            output["aligned_valid_reads"] = 1
            
            if output["flag"] == flag_list[0]: #first read in pair oriented 5'to 3' (positive)
                output["orientation"]  = "+"
                output["border"] = output["sequence"][:2] 

            elif output["flag"] == flag_list[1]: #first read in pair oriented 3'to 5' (negative)
                output["orientation"] = "-"
                output["border"] = output["sequence"][::-1][:2] #needs to be reversed to make sure the start position is always the same

                #for CIGAR
                matches = findall(r'(\d+)([A-Z]{1})', output["cigar"])
                output["clipped"] = 0
                for match in matches:
                    if match[1] == "S":
                        output["clipped"]=int(match[0])
                        break #only consideres the first one at the start

                output["local"]=str(int(output["local"])+len(output["sequence"])-output["clipped"]-1) # -1 to offsset bowtie alignement

    return output