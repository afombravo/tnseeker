import multiprocessing
from datetime import datetime
from colorama import Fore
import csv
import subprocess
import os
import glob

def cpu():
    c = multiprocessing.cpu_count()
    if c >= 2:
        c -= 1
    return c

def colourful_errors(warning_type:str, error: str):
    warning_colour = Fore.GREEN
    if warning_type == "FATAL":
        warning_colour = Fore.RED
    elif warning_type == "WARNING":
        warning_colour = Fore.YELLOW

    print(f"{Fore.BLUE} {datetime.now().strftime('%c')}{Fore.RESET} [{warning_colour}{warning_type}{Fore.RESET}] {error}")

def csv_writer(output_file_path,output_file):
    with open(output_file_path, "w", newline='') as output:
        writer = csv.writer(output)
        writer.writerows(output_file)

def subprocess_cmd(command: list):
    try:
        return subprocess.check_output(command)
    except subprocess.CalledProcessError as e:
        return e.output.decode()
    
def file_finder(folder: str, ext: list, search_term: str = None):
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
    variables = {}
    with open(var_file) as current:
        for line in current:
            key,value = line[:-1].split(" : ")
            variables[key] = value
    return variables

def adjust_spines(ax, spines,x,y): #offset spines
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 2.5))  # outward by 10 points
            if (loc == 'left') or (loc == 'right'):
                spine.set_bounds(y)
            else:
                spine.set_bounds(x)
        else:
            spine.set_color('none')  # don't draw spine

def gb_parser(file_line):
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