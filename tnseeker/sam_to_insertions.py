import csv
import numpy as np
import os, glob
from regex import findall
from matplotlib import pyplot as plt
import argparse
from Bio import SeqIO
import multiprocessing
import pkg_resources
import subprocess
import datetime
from colorama import Fore

""" This script processes and analyzes sequencing data from a SAM (Sequence Alignment/Map) file. 
    The main purpose is to extract information about transposon insertions in a given genome 
    and generate statistics, including the read histogram, for further analysis.
    
    Here's a high-level overview of the code:

        Define the main function, which takes a list of command line arguments 
        and processes them. It sets up the input parameters, such as the folder 
        path, name folder, paired-ended or single-ended reads, read threshold, 
        barcode, and map quality threshold. Then it calls the path_finder and 
        extractor functions.
        
        The path_finder function finds all SAM files in the given folder and 
        returns a list of their file paths.
        
        The extractor function processes each SAM file, filtering the aligned 
        reads based on the specified parameters (flags, map quality threshold, etc.). 
        It then builds a dictionary (insertion_count) that stores information 
        about unique Tn5 insertions, including their location, orientation, 
        border sequence, read count, and average mapping quality. If barcodes 
        are used, it will also store barcode-related information.
        
        The plotter function generates a histogram of read distribution based 
        on the unique Tn5 insertions and saves it as a PNG image.
        
        The dictionary_parser and insert_parser functions write the transposon insertion 
        information into CSV files. The former generates a CSV file containing all 
        unique insertions, while the latter generates a CSV file containing 
        barcode-specific information (if barcodes are used).
        
        The main function writes a summary of the analysis results, including 
        total aligned reads, quality passed reads, unique insertions, and read 
        distribution statistics, to a text file.
        
        Finally, the script runs the main function with the command line arguments provided.
    
    In summary, this script processes sequencing data, extracts transposon insertion 
    information, generates a read histogram, and outputs the results in text and 
    CSV files for further analysis."""
    
def main(argv):
    folder_path = argv[0]
    name_folder = argv[1]
    paired_ended = argv[2]
    read_threshold = argv[3] == "True"
    read_cut = int(argv[4]) if read_threshold else 0
    barcode = argv[5] == "True"
    map_quality_threshold = int(argv[6])
    annotation_file = argv[7]
    ir_size_cutoff = int(argv[8])
    
    pathing = path_finder(folder_path)
    extractor(name_folder, folder_path, pathing, paired_ended,barcode,\
              read_threshold,read_cut,annotation_file,ir_size_cutoff,map_quality_threshold)

def path_finder(folder_path): 
    filenames = []
    for filename in glob.glob(os.path.join(folder_path, '*.sam')):
        filenames.append(filename)
    return filenames

def getNextFilePath(output_folder, name, output_file):

    output_file_path = os.path.join(output_folder, name + ".csv")
    
    with open(output_file_path, "w", newline='') as output:
        writer = csv.writer(output)
        writer.writerows(output_file)

def adjust_spines(ax, spines,x,y): #offset spines
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))  # outward by 10 points
            if loc == 'left':
                spine.set_bounds(y)
            else:
                spine.set_bounds(x)
        else:
            spine.set_color('none')  # don't draw spine

def plotter(insertion_count, naming, output_folder):
    
    reads = []
    for key in insertion_count:
        reads.append(insertion_count[key].count)

    log_reads = np.log10(reads)
    
    median = np.median(reads)
    average = int(np.average(reads))
    std = int(np.std(reads))
    
    statstics = f"Read Distribution for {naming}: Median: {median}; Average: {average}; Std: {std}\n"
    
    fig, ax1 = plt.subplots()  

    plt.hist(log_reads, weights=np.ones(len(log_reads)) / len(log_reads), bins = 30, color = "orange")

    plt.title('Read Histogram')
    plt.xlabel("Reads per insertion (log10)")
    plt.ylabel("Tn5 Insertion Frequency")

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(True)
    ax1.spines['left'].set_visible(True)
    
    ax1.legend([naming], loc="upper right",prop={'size': 6.5})
    
    plt.savefig(output_folder + f"/read_distribution_{naming}.png", dpi=300)
    plt.close()
    
    return statstics

class Insertion():
    
    def __init__(self, contig=None, local=None, border=None, orientation=None,\
                 count=None,barcode=None,name=None,product=None,gene_orient=None,
                 relative_gene_pos=None,mapQ=0,read_id=None):
        
        self.contig = contig
        self.local = local
        self.orientation = orientation
        self.seq = border
        self.count = count
        self.mapQ = mapQ
        self.barcode = barcode or {}
        self.name = name
        self.product = product
        self.gene_orient = gene_orient
        self.relative_gene_pos = relative_gene_pos
        self.read_id = read_id

def cpu():
    c = multiprocessing.cpu_count()
    if c >= 2:
        c -= 1
    pool = multiprocessing.Pool(processes = c)
    return pool, c

def subprocess_cmd(command):
    try:
        return subprocess.check_output(command)
    except subprocess.CalledProcessError as e:
        return e.output.decode()

def extractor(name_folder, folder_path, pathing, paired_ended,barcode,\
              read_threshold,read_cut,annotation_file,ir_size_cutoff,map_quality_threshold = 42):
    
    def barcode_finder():
        read = ""
        barcode = {}
        bar=None
        with open(f"{folder_path}/barcodes_1.txt") as current:
            for line in current:
                if "@" not in line:
                    if len(line[:-1]) != 0:
                        bar = line[:-1]
                    else:
                        bar = None
                        
                if (bar is not None) & ("@" in line):
                    read = line[:-1].split(" ")[0][1:]
                    barcode[read] = bar
                    bar = None

        return barcode

    aligned_reads, aligned_valid_reads = 0, 0
    insertion_count = {}
    
    flag_list = [0, 16]
    if paired_ended=="PE":
        flag_list = [83, 99] #[16] for single ended data #99 and 83 means that the read is the first in pair (only paired ended reads are considered as valid)
    
    file = pathing[0]
    if barcode:
        print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.GREEN}INFO{Fore.RESET}] Linking barcodes to insertion positions")
        file = f"{folder_path}/barcoded_align.sam"
        subprocess_cmd([pkg_resources.resource_filename(__name__, 'data/barcode2sam.sh'),
                       f"{folder_path}/alignment.sam",
                       f"{folder_path}/barcodes_1.txt",
                       file,
                       folder_path])
    
    print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.GREEN}INFO{Fore.RESET}] Parsing Bowtie alignments into an insertion matrix")
    with open(file) as current:
        for line in current:
            sam = line.split('\t')
            if (sam[0][0] != "@") and (sam[2] != '*'): #ignores headers and unaligned contigs
                local = sam[3]
                sequence = sam[9] 
                contig = sam[2]
                flag = int(sam[1])
                cigar = sam[5]
                map_quality = float(sam[4])
                border, orientation = "",""
                aligned_reads += 1
                multi = "XS:i:" in sam #multiple alignemnts

                if (flag in flag_list) & (multi==False) & (map_quality >= map_quality_threshold): #only returns aligned reads witht he proper flag score
                    
                    aligned_valid_reads += 1
                    
                    if flag == flag_list[0]: #first read in pair oriented 5'to 3' (positive)
                        orientation = "+"
                        border = sequence[:2] 

                    elif flag == flag_list[1]: #first read in pair oriented 3'to 5' (negative)
                        orientation = "-"
                        border = sequence[::-1][:2] #needs to be reversed to make sure the start position is always the same

                        #for CIGAR
                        matches = findall(r'(\d+)([A-Z]{1})', cigar)
                        clipped = 0
                        for match in matches:
                            if match[1] == "S":
                                clipped=int(match[0])
                                break #only consideres the first one at the start

                        local=str(int(local)+len(sequence)-clipped-1) # -1 to offsset bowtie alignement

                    key = (contig, local, orientation)
                    if key not in insertion_count: 
                        insertion_count[key] = Insertion(contig=key[0], 
                                                         local=key[1], 
                                                         orientation=key[2], 
                                                         count=1, 
                                                         border=border,
                                                         mapQ=map_quality)
                    else: 
                        insertion_count[key].count += 1
                        insertion_count[key].mapQ += map_quality
                    
                    if barcode:
                        bar = None
                        if "BC:Z:" in sam[-1]:
                            bar = sam[-1][:-1].split(":")[2]
                        if bar != None:
                            if bar in insertion_count[key].barcode:
                                insertion_count[key].barcode[bar] += 1
                            else:
                                insertion_count[key].barcode[bar] = 1

    for key in insertion_count:
        insertion_count[key].mapQ = insertion_count[key].mapQ / insertion_count[key].count    

    pool,cpus = cpu()
    len_insertion_count_divider = int(len(insertion_count) / cpus)
    batch_goals = {}
    for i in range(cpus):
        batch_goals[i] = set()
        for j,key in enumerate(insertion_count):
            if j >= len_insertion_count_divider * i:
                if i != cpus-1:
                    if len(batch_goals[i]) <= len_insertion_count_divider:
                        batch_goals[i].add(key)
                    else:
                        break
                else:
                    batch_goals[i].add(key)

    result_objs = []
    for batch in batch_goals:

        insertion_count_filtered = {}
        for key in insertion_count:
            if key in batch_goals[batch]:
                insertion_count_filtered[key] = insertion_count[key]

        result=pool.apply_async(annotation_processer, 
                            args=((insertion_count_filtered, 
                                   read_threshold,
                                   read_cut,
                                   barcode,
                                   annotation_file,
                                   ir_size_cutoff,
                                   name_folder,
                                   folder_path)))
    
        result_objs.append(result)
    pool.close()
    pool.join()
        
    result = [result.get() for result in result_objs]
    
    final_compiler = {}
    barcoded_insertions_final = []
    insertions_final = []

    for entry in result:
        insertion,barcoded,insert = entry
        final_compiler = {**final_compiler, **insertion}

        if barcode:
            for barcode in barcoded:
                barcoded_insertions_final.append(barcode)
        
        for insertion in insert:
            insertions_final.append(insertion)
    
    if barcode:
        annotate_barcodes_writer(barcoded_insertions_final,insertions_final,name_folder,folder_path)
        
    dictionary_parser(final_compiler,barcode,folder_path,name_folder)
        
    q = plotter(insertion_count, f"Unique insertions_{name_folder}", folder_path)

    reads = f" Total Aligned Reads: {aligned_reads}\nTotal Quality Passed Reads: {aligned_valid_reads}\nFiltered Vs. Raw Read % ratio: {round(aligned_valid_reads/aligned_reads*100,2)}%\n"
    e = " Number of total unique insertions: {}\n".format(len(insertion_count))
    
    print(f"\n{Fore.YELLOW} -- Library statistics -- {Fore.RESET}\n")
    print(f"{Fore.GREEN} Total aligned reads: {Fore.RESET}{aligned_reads}")
    print(f"{Fore.GREEN} Total quality passed reads: {Fore.RESET}{aligned_valid_reads}")
    print(f"{Fore.GREEN} Filtered Vs. Raw Read % ratio: {Fore.RESET}{round(aligned_valid_reads/aligned_reads*100,2)}%")
    print(f"{Fore.GREEN} Number of total unique insertions: {Fore.RESET}{len(insertion_count)}")
    print(f"\n{Fore.YELLOW} ---- {Fore.RESET}\n")
    
    with open("{}/library_stats_{}.txt".format(folder_path,name_folder), "w+") as current:
        current.write(reads+q+e)

def annotation_processer(insertion_count_filtered,read_threshold,read_cut,
                         barcode,annotation_file,ir_size_cutoff,name_folder,folder_path):

    if read_threshold:
        insertion_count_filtered=dict_filter(insertion_count_filtered,read_cut,barcode)
    
    if (annotation_file.endswith(".gb")) or (annotation_file.endswith(".gbk")):
        insertion_count_filtered,genes,contigs = gene_parser_genbank(annotation_file,insertion_count_filtered)
        
    elif annotation_file.endswith(".gff"):
        insertion_count_filtered,genes,contigs = gene_parser_gff(annotation_file,insertion_count_filtered)
        
    insertion_count_filtered = inter_gene_annotater(annotation_file,insertion_count_filtered,ir_size_cutoff,genes,contigs)
    
    barcoded_insertions,insertions = [],[]
    if barcode:
        barcoded_insertions,insertions = insert_parser(insertion_count_filtered,name_folder,folder_path,barcode)

    return insertion_count_filtered,barcoded_insertions,insertions

def dict_filter(dictionary,read_cut):
    for key in list(dictionary):
        if dictionary[key].count < read_cut:
            del dictionary[key]
    return dictionary

def insert_parser(insertion_count,name_folder,folder_path,barcode):    
    insertions,barcoded_insertions = [],[]

    for key in insertion_count: 

        contig = [insertion_count[key].contig]
        local = [insertion_count[key].local]
        orientation = [insertion_count[key].orientation]
        count = [insertion_count[key].count]
        mapq = [insertion_count[key].mapQ]
        gene_name = [insertion_count[key].name]
        gene_product = [insertion_count[key].product]
        gene_orientation = [insertion_count[key].gene_orient]
        relative_gene_pos = [insertion_count[key].relative_gene_pos]
        
        barcodes,reads = '',0
        for bar,read in insertion_count[key].barcode.items():
            barcodes += f'{bar}:{read};'
            reads += read
            
            ## for individual barcoded insertions
            barcoded_insertions.append([bar] + [read] + contig + local +\
                                       orientation + count + mapq + gene_name + \
                                       gene_product + gene_orientation + relative_gene_pos)
    

        insertions.append(contig + local + orientation + count + mapq + \
                          gene_name + gene_product + gene_orientation + relative_gene_pos + \
                          [len(insertion_count[key].barcode)] + [reads] + [barcodes])

    return barcoded_insertions,insertions

def annotate_barcodes_writer(barcoded_insertions,insertions,name_folder,folder_path):
    
    insertions.insert(0, ["#Contig"] + ["position"] + ["Orientation"] + ["Total Reads"] + \
                      ["Average MapQ"] + ["Gene Name"] + ["Gene Product"] + ["Gene Orientation"] + \
                      ["Relative Position in Gene (0-1)"] + ["Number of different barcodes in coordinate"] + \
                      ["Total barcode Reads"] + ["Barcodes (barcode:read)"])
    
    name = f"barcoded_insertions_{name_folder}.csv"
    output_file_path = os.path.join(folder_path, name)
    with open(output_file_path, "w", newline='') as output:
        writer = csv.writer(output)
        writer.writerows(insertions)
    
    ############
    
    barcoded_insertions.insert(0, ["#Barcode"] + ["Barcode Reads"] +\
                               ["Contig"] + ["position"] + ["Orientation"] + ["Total Reads in position"] + \
                              ["Average MapQ"] + ["Gene Name"] + ["Gene Product"] + ["Gene Orientation"] + \
                              ["Relative Position in Gene (0-1)"])
        
    name = f"annotated_barcodes_{name_folder}.csv"
    output_file_path = os.path.join(folder_path, name)
    with open(output_file_path, "w", newline='') as output:
        writer = csv.writer(output)
        writer.writerows(barcoded_insertions)

def inter_gene_annotater(annotation_file,insertion_count,ir_size_cutoff,genes,contigs):

    ir_annotation = {}
    count = 0
    for i,gene in enumerate(genes[:-1]):
        contig = gene[-1]
        if contig == genes[i+1][-1]: #same contigs
            gene_down_start_border = ir_size_cutoff + gene[1]
            gene_up_start_border = genes[i+1][0] - ir_size_cutoff
            domain_size = gene_up_start_border - gene_down_start_border
            if domain_size >= 1:
                count += 1
                ir_annotation[f'IR_{count}_{gene[3]}_{gene[2]}_UNTIL_{genes[i+1][3]}_{gene[2]}'] = (gene[1],genes[i+1][0],domain_size,contig)

        if contig != genes[i+1][-1]:

            circle_closer = gene[1] + ir_size_cutoff 
            domain_size = contigs[contig] - circle_closer
            if domain_size >= 1:
                count += 1
                ir_name = f'IR_{count}_{gene[3]}_{gene[2]}_contig_{contig}_-end'
                if ir_name not in ir_annotation:
                    ir_annotation[ir_name] = (genes[-1][1],contigs[contig],domain_size,contig)

    circle_closer = genes[-1][1] + ir_size_cutoff 
    domain_size = contigs[contig] - circle_closer
    if domain_size >= 1:
        count += 1
        ir_name = f'IR_{count}_{genes[-1][3]}_{genes[-1][2]}_contig_{contig}_-end'
        if ir_name not in ir_annotation:
            ir_annotation[ir_name] = (genes[-1][1],contigs[contig],domain_size,contig)

    domain_size = genes[0][0] - ir_size_cutoff
    if domain_size  >= 1:
        count += 1
        ir_name = f'IR_{count}_contig_{contig}_-start_{genes[0][3]}_{genes[0][2]}'
        if ir_name not in ir_annotation:
           ir_annotation[ir_name] = (0,genes[0][0],domain_size,contig)
    
    for ir in ir_annotation:
        for key in insertion_count:
            if insertion_count[key].name is None:
                if insertion_count[key].contig == ir_annotation[ir][3]:
                    if (int(insertion_count[key].local) >= ir_annotation[ir][0]) & (int(insertion_count[key].local) <= ir_annotation[ir][1]):
                        insertion_count[key].name = ir
                        insertion_count[key].relative_gene_pos = (int(insertion_count[key].local) - ir_annotation[ir][0]) / ir_annotation[ir][2]
             
    return insertion_count

def gene_parser_genbank(annotation_file,insertion_count):
    
    ''' The gene_info_parser_genbank function takes a genbank file as input and 
    extracts gene information, storing it in a dictionary with Gene class 
    instances as values. It parses the file using the SeqIO module, 
    retrieving attributes such as start, end, orientation, identity, 
    and product for each gene.''' 
    
    contigs = {}
    genes = []
    for rec in SeqIO.parse(annotation_file, "gb"):
        for feature in rec.features:
            if feature.type != 'source':
                start = feature.location.start
                end = feature.location.end
                domain_size = end - start
                
                orientation = feature.location.strand
                
                try:
                    identity = feature.qualifiers['locus_tag'][0]
                except KeyError:
                    identity = None

                if orientation == 1:
                    orientation = "+"
                else:
                    orientation = "-"
                
                try:
                    if 'product' in feature.qualifiers:
                        product = feature.qualifiers['product'][0]
                    else:
                        product = feature.qualifiers['note'][0]
                except KeyError:
                    product = None
                    
                for key, val in feature.qualifiers.items():   
                    if "pseudogene" in key:
                        gene = identity
                        break #avoids continuing the iteration and passing to another key, which would make "gene" assume another value
                    elif "gene" in key:
                        gene = feature.qualifiers['gene'][0]
                        break #avoids continuing the iteration and passing to another key, which would make "gene" assume another value
                    else:
                        gene = identity

                genes.append((start,end,orientation,gene,rec.id))

                for key in insertion_count:
                    if insertion_count[key].contig == rec.id:
                        if (int(insertion_count[key].local) >= start) & (int(insertion_count[key].local) <= end):
                            insertion_count[key].name = gene
                            insertion_count[key].product = product
                            insertion_count[key].gene_orient = orientation
                            insertion_count[key].relative_gene_pos = (int(insertion_count[key].local) - start) / domain_size

        contigs[rec.id] = len(rec.seq)
        genes = list(dict.fromkeys(genes))
        genes.sort(key=lambda x: (x[-1], x[0])) #sort by start position of the gene and contig
    return insertion_count,genes,contigs

def gene_parser_gff(annotation_file,insertion_count):
    
    contigs = {}
    genes = []
    with open(annotation_file) as current:
        for line in current:
            GB = line.split('\t') #len(GB)
            
            if "##FASTA" in GB[0]:
                break
            
            if "#" not in GB[0][:3]: #ignores headers

                start = int(GB[3])
                end = int(GB[4])
                domain_size = end - start
                
                features = GB[8].split(";") #gene annotation file
                feature = {}
                for entry in features:
                    entry = entry.split("=")
                    if len(entry) == 2:
                        feature[entry[0]] = entry[1].replace("\n","")
                if "gene" in feature:
                    gene=feature["gene"]
                if "Name" in feature:
                    gene=feature["Name"]
                else:
                    gene=feature["ID"]
                    
                if 'product' not in feature:
                    feature['product'] = None
                    
                if gene == ".":
                    gene = f"{GB[2]}_{start}_{end}"
                
                contig = GB[0]
                orientation = GB[6] #orientation of the gene
                
                for key in insertion_count:
                    if insertion_count[key].contig == contig:
                        if (int(insertion_count[key].local) >= start) & (int(insertion_count[key].local) <= end):
                            insertion_count[key].name = gene
                            insertion_count[key].product =  feature['product']
                            insertion_count[key].gene_orient = orientation
                            insertion_count[key].relative_gene_pos = (int(insertion_count[key].local) - start) / domain_size
                
                key = (start,end,orientation,gene,contig)
                genes.append(key)

            if "sequence-region" in GB[0]:
                GB = GB[0].split(" ")
                contigs[GB[-3]] = int(GB[-1][:-1])
                
        genes = list(dict.fromkeys(genes))
        genes.sort(key=lambda x: (x[-1], x[0])) #sort by start position of the gene and contig
        
        if len(genes) == 0:
            print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.YELLOW}WARNING{Fore.RESET}] Watch out, no genomic features were loaded. The gff file is not being parsed correctly, or was not loaded.")
        
    return insertion_count,genes,contigs

def dictionary_parser(dictionary,barcode,folder_path,name_folder):
    
    insertions = []

    for key in dictionary:

        contig = [dictionary[key].contig]
        local = [dictionary[key].local]
        orientation = [dictionary[key].orientation]
        border = [dictionary[key].seq]
        count = [dictionary[key].count]
        gene_name = [dictionary[key].name]
        gene_product = [dictionary[key].product]
        gene_orientation = [dictionary[key].gene_orient]
        relative_gene_pos = [dictionary[key].relative_gene_pos]
        mapQ = [dictionary[key].mapQ]
        
        insertions.append(contig + local + orientation + border + count + mapQ +\
                          gene_name + gene_product + gene_orientation + relative_gene_pos)

    insertions.insert(0, ["#Contig"] + ["position"] + ["Orientation"] + \
                      ["Transposon Border Sequence"] + ["Read Counts"] + \
                      ["Average mapQ across reads"] + ["Gene Name"] + ["Gene Product"] + \
                      ["Gene Orientation"] + ["Relative Position in Gene (0-1)"])

    getNextFilePath(folder_path, "all_insertions_"+name_folder, insertions) #all the unique insertions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process SAM aligned files and extract relevant information.")
    parser.add_argument("folder_path", help="Path to the folder containing SAM files.")
    parser.add_argument("name_folder", help="Name of the folder.")
    parser.add_argument("paired_ended", help="Type of data: 'PE' for paired-end or 'SE' for single-end.")
    parser.add_argument("read_threshold", type=bool, help="Apply read threshold (True/False).")
    parser.add_argument("read_cut", type=int, help="Read cut value.")
    parser.add_argument("barcode", type=bool, help="Use barcodes (True/False).")
    parser.add_argument("map_quality_threshold", type=int, help="Map quality threshold.")
    parser.add_argument("gb_annotation_file", help="Needs to be a standard .gb file")
    parser.add_argument("ir_size_cutoff", type=int, help="The number of bp up and down stream of any gene to be considered an intergenic region")

    args = parser.parse_args()
    
    main([args.folder_path, 
          args.name_folder, 
          args.paired_ended, 
          args.read_threshold, 
          args.read_cut, 
          args.barcode, 
          args.map_quality_threshold,
          args.gb_annotation_file,
          args.ir_size_cutoff])
