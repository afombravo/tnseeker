import numpy as np
import os
from tnseeker.extras.helper_functions import colourful_errors,csv_writer,file_finder,variables_parser,gb_parser,gff_parser,inter_gene_annotater,bowtie2parser
from matplotlib import pyplot as plt
from Bio import SeqIO
import multiprocessing
from colorama import Fore
import sys

""" This script processes and analyzes sequencing data from a SAM (Sequence Alignment/Map) file. 
    The main purpose is to extract information about transposon insertions in a given genome 
    and generate statistics, including the read histogram, for further analysis.
"""

def main(argv):
    global variables
    variables = variables_parser(argv)
    variables['files'] = file_finder(variables['directory'], ['*.sam'])
    variables["read_threshold"] = variables["read_threshold"] == "True"
    variables["read_value"] = int(variables["read_value"]) if variables["read_threshold"] else 0
    variables["barcode"] = variables["barcode"] == "True"
    variables["MAPQ"] = int(variables["MAPQ"])
    variables['intergenic_size_cutoff'] = int(variables['intergenic_size_cutoff'])
    variables['cpus'] = int(variables['cpus'])

    extractor()

def plotter(insertion_count, naming, output_folder):
    
    reads = []
    for key in insertion_count:
        reads.append(insertion_count[key].count)

    log_reads = np.log10(reads)
    
    median = np.median(reads)
    average = int(np.average(reads))
    std = int(np.std(reads))
    
    statstics = f"Read distribution for {naming}:\n  Median: {median}\n  Average: {average}\n  Std: {std}\n"
    
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

def extractor():

    aligned_reads, aligned_valid_reads = 0, 0
    insertion_count = {}
    
    flag_list = [0, 16]
    if variables['seq_type']=="PE":
        flag_list = [83, 99] 
        
    sam_file = variables['files'][0]

    colourful_errors("INFO",
        "Parsing Bowtie alignments into an insertion matrix.")
    
    with open(sam_file) as current:
        for line in current:
            sam = line.split('\t')
            sam_output = bowtie2parser(sam,variables["MAPQ"],flag_list)
            
            if "aligned_reads" in sam_output:
                aligned_reads += 1

            if "aligned_valid_reads" in sam_output:
                aligned_valid_reads += 1
                key = (sam_output["contig"], sam_output["local"], sam_output["orientation"])
                if key not in insertion_count: 
                    insertion_count[key] = Insertion(contig=sam_output["contig"], 
                                                        local=sam_output["local"], 
                                                        orientation=sam_output["orientation"], 
                                                        count=1, 
                                                        border=sam_output["border"],
                                                        mapQ=sam_output["map_quality"])
                else: 
                    insertion_count[key].count += 1
                    insertion_count[key].mapQ += sam_output["map_quality"]
                
                if variables['barcode']:
                    bar = None
                    if ":BC:" in sam[0]:
                        bar = sam[0].split(":BC:" )[1]
                    if bar != None:
                        if bar in insertion_count[key].barcode:
                            insertion_count[key].barcode[bar] += 1
                        else:
                            insertion_count[key].barcode[bar] = 1

    for key in insertion_count:
        insertion_count[key].mapQ = insertion_count[key].mapQ / insertion_count[key].count    

    len_insertion_count_divider = int(len(insertion_count) / variables['cpus'])
    batch_goals = {}
    for i in range(variables['cpus']):
        batch_goals[i] = set()
        for j,key in enumerate(insertion_count):
            if j >= len_insertion_count_divider * i:
                if i != variables['cpus']-1:
                    if len(batch_goals[i]) <= len_insertion_count_divider:
                        batch_goals[i].add(key)
                    else:
                        break
                else:
                    batch_goals[i].add(key)

    colourful_errors("INFO","Annotating insertions.")

    annotations_dict,contigs = annotation_ETL()

    if not annotations_dict:
        colourful_errors("WARNING","Watch out, no genomic features were found. The annotation file is not being parsed correctly, or was not loaded - check paths.")

    result_objs = []
    pool = multiprocessing.Pool(processes = variables['cpus'])
    for batch in batch_goals:

        insertion_count_filtered = {}
        for key in insertion_count:
            if key in batch_goals[batch]:
                insertion_count_filtered[key] = insertion_count[key]

        result=pool.apply_async(annotation_processer, 
                            args=(annotations_dict,insertion_count_filtered))
    
        result_objs.append(result)
    pool.close()
    pool.join()
        
    result = [result.get() for result in result_objs]
    
    final_compiled_insertions = {}
    for insertion in result:
        final_compiled_insertions = {**final_compiled_insertions, **insertion}

    intergenic_regions = inter_gene_annotater(annotations_dict, contigs, variables['intergenic_size_cutoff'])

    for ir in intergenic_regions:
        for key in final_compiled_insertions:
            if final_compiled_insertions[key].name is None:
                if final_compiled_insertions[key].contig == intergenic_regions[ir][3]:
                    if intergenic_regions[ir][0] <= int(final_compiled_insertions[key].local) <= intergenic_regions[ir][1]:
                        final_compiled_insertions[key].name = ir
                        final_compiled_insertions[key].relative_gene_pos = (int(final_compiled_insertions[key].local) - intergenic_regions[ir][0]) / intergenic_regions[ir][2]

    if variables['barcode']:
        barcoded_insertions_final = []
        insertions_barcoded_final = []
        barcoded2insertions,insertions2barcode = barcode_demultiplexer(final_compiled_insertions)

        for barcode in barcoded2insertions:
            barcoded_insertions_final.append(barcode)
        
        for insertion in insertions2barcode:
            insertions_barcoded_final.append(insertion)

        barcoded_insertions_final.sort(key=lambda x: (x[2], int(x[3])))
        barcoded_insertions_final.insert(0,["#Barcode"] + ["Barcode Reads"] +\
                                            ["Contig"] + ["position"] + ["Orientation"] + ["Total Reads in position"] + \
                                            ["Average MapQ"] + ["Gene Name"] + ["Gene Product"] + ["Gene Orientation"] + \
                                            ["Relative Position in Gene (0-1)"])
        
        name = f"annotated_barcodes_{variables['strain']}.csv"
        output_file_path = os.path.join(variables['directory'], name)
        csv_writer(output_file_path,barcoded_insertions_final)

        insertions_barcoded_final.sort(key=lambda x: (x[0], int(x[1])))
        insertions_barcoded_final.insert(0,["#Contig"] + ["position"] + ["Orientation"] + ["Total Reads"] + \
                                            ["Average MapQ"] + ["Gene Name"] + ["Gene Product"] + ["Gene Orientation"] + \
                                            ["Relative Position in Gene (0-1)"] + ["Number of different barcodes in coordinate"] + \
                                            ["Total barcode Reads"] + ["Barcodes (barcode:read)"])
        name = f"barcoded_insertions_{variables['strain']}.csv"
        output_file_path = os.path.join(variables['directory'], name)
        csv_writer(output_file_path,insertions_barcoded_final)

    insertion_demultiplexer(final_compiled_insertions)
        
    q = plotter(insertion_count, f"Unique insertions_{variables['strain']}", variables['directory'])
    reads = f"""\n    ####  \nALIGNMENT INFO\n    ####  \n\nTotal detected aligned reads: {aligned_reads}\nTotal quality passed reads: {aligned_valid_reads}\nQuality passing Vs. total aligned reads % ratio: {round(aligned_valid_reads/aligned_reads*100,2)}%\n"""
    total_insertions_count = f"Number of unique insertions in the library: {len(insertion_count)}\n"
    
    print(f"\n{Fore.YELLOW} -- Library statistics -- {Fore.RESET}\n")
    print(f"{Fore.GREEN} Total aligned reads: {Fore.RESET}{aligned_reads}")
    print(f"{Fore.GREEN} Total quality passed reads: {Fore.RESET}{aligned_valid_reads}")
    print(f"{Fore.GREEN} Quality passing Vs. total aligned reads % ratio: {Fore.RESET}{round(aligned_valid_reads/aligned_reads*100,2)}%")
    print(f"{Fore.GREEN} Number of total unique insertions: {Fore.RESET}{len(insertion_count)}")
    print(f"\n{Fore.YELLOW} ---- {Fore.RESET}\n")
    
    with open("{}/library_stats_{}.log".format(variables['directory'],variables['strain']), "w+") as current:
        current.write(reads+q+total_insertions_count)

def annotation_processer(annotations_dict, insertion_count_filtered):

    if variables['read_threshold']:
        insertion_count_filtered=dict_filter(insertion_count_filtered,
                                             variables['read_value'])

    return insertion_annotator(annotations_dict, insertion_count_filtered)

def annotation_ETL():

    if (variables['annotation_file'].endswith(".gb")) or (variables['annotation_file'].endswith(".gbk")):
        return gene_parser_genbank()
        
    elif variables['annotation_file'].endswith(".gff"):
        return gene_parser_gff()

def dict_filter(dictionary,read_cut):
    for key in list(dictionary):
        if dictionary[key].count <= read_cut:
            del dictionary[key]
    return dictionary

def barcode_demultiplexer(insertion_count):    
    insertions2barcode,barcoded2insertions = [],[]

    for key in insertion_count: 
        barcodes,reads = '',0
        for bar,read in insertion_count[key].barcode.items():
            barcodes += f'{bar}:{read};'
            reads += read
            barcoded2insertions.append([bar] + \
                                       [read] + \
                                       [insertion_count[key].contig] + \
                                       [insertion_count[key].local] + \
                                       [insertion_count[key].orientation] + \
                                       [insertion_count[key].count] + \
                                       [insertion_count[key].mapQ] + \
                                       [insertion_count[key].name] + \
                                       [insertion_count[key].product] + \
                                       [insertion_count[key].gene_orient] + \
                                       [insertion_count[key].relative_gene_pos])
    
        insertions2barcode.append([insertion_count[key].contig] + \
                          [insertion_count[key].local] + \
                          [insertion_count[key].orientation] + \
                          [insertion_count[key].count] + \
                          [insertion_count[key].mapQ] + \
                          [insertion_count[key].name] + \
                          [insertion_count[key].product] + \
                          [insertion_count[key].gene_orient] + \
                          [insertion_count[key].relative_gene_pos] + \
                          [len(insertion_count[key].barcode)] + \
                          [reads] + \
                          [barcodes])

    return barcoded2insertions,insertions2barcode

def insertion_annotator(annotations_dict, insertion_count):

    for gene_info in annotations_dict.values():
        start, end, orientation = gene_info["start"], gene_info["end"], gene_info["orientation"]

        for insertion in insertion_count.values():
            if insertion.contig == gene_info["contig"]:
                local_pos = int(insertion.local)
                if start <= local_pos <= end:
                    insertion.name = gene_info["gene"]
                    insertion.product = gene_info["product"]
                    insertion.gene_orient = orientation

                    relative_pos = (local_pos - start) / gene_info["domain_size"]
                    insertion.relative_gene_pos = 1 - relative_pos if orientation == "-" else relative_pos

    return insertion_count

def gene_parser_genbank():

    contigs = {}
    genes = {}
    for rec in SeqIO.parse(variables['annotation_file'], "gb"):
        for feature in rec.features:
            if feature.type != 'source':
                gene_info = gb_parser(feature)

                if gene_info:
                    gene_info["contig"] = rec.id
                    gene_key = ((gene_info["start"], gene_info["end"], gene_info["orientation"], gene_info["gene"], rec.id))
                    genes[gene_key] = gene_info

        contigs[rec.id] = len(rec.seq)

    if not contigs:
        colourful_errors("FATAL","No headers with contig names and sizes were found in the genbank file. Please check the file.")
        exit()

    return genes,contigs

def gene_parser_gff():
    
    contigs = {}
    genes = {}
    with open(variables['annotation_file']) as current:
        for line in current:
            line = line.split('\t')
            gene_info = gff_parser(line)

            if gene_info:
                gene_key = ((gene_info["start"], gene_info["end"], gene_info["orientation"], gene_info["gene"], gene_info["contig"]))
                genes[gene_key] = gene_info

            if "sequence-region" in line[0]:
                line = line[0].split(" ")
                contigs[line[-3]] = int(line[-1][:-1])

    if not contigs:
        colourful_errors("FATAL","No headers with contig names and sizes were found in the GFF file. Please check the headers of the example gff file at github for reference.")
        exit()
    
    return genes,contigs

def insertion_demultiplexer(dictionary):
    
    insertions = []
    for key in dictionary:
        insertions.append([dictionary[key].contig] + \
                          [dictionary[key].local] + \
                          [dictionary[key].orientation] + \
                          [dictionary[key].seq] + \
                          [dictionary[key].count] + \
                          [dictionary[key].mapQ] + \
                          [dictionary[key].name] + \
                          [dictionary[key].product] + \
                          [dictionary[key].gene_orient] + \
                          [dictionary[key].relative_gene_pos])

    insertions.sort(key=lambda x: (x[0], int(x[1])))

    bedgraph_plus,bedgraph_minus = "",""
    for key in insertions:
        if key[2] == "+":
            bedgraph_plus += f"{key[0]}\t{key[1]}\t{int(key[1])+1}\t{key[4]}\n"
        else:
            bedgraph_minus += f"{key[0]}\t{key[1]}\t{int(key[1])+1}\t{key[4]}\n"

    insertions.insert(0,["#Contig"] + ["position"] + ["Orientation"] + \
                        ["Transposon chromosome Border Sequence"] + ["Read Counts"] + \
                        ["Average mapQ across reads"] + ["Gene Name"] + ["Gene Product"] + \
                        ["Gene Orientation"] + ["Relative Position in Gene (0-1)"])
    
    output_file_path = os.path.join(variables['directory'], f"all_insertions_{variables['strain']}.csv") #all the unique insertions
    csv_writer(output_file_path,insertions)

    #write the bedgraph file
    with open(os.path.join(variables['directory'], f"all_insertions_{variables['strain']}_bedgraph_plus.bedgraph"), "w+") as current:
        current.write(bedgraph_plus)
    with open(os.path.join(variables['directory'], f"all_insertions_{variables['strain']}_bedgraph_minus.bedgraph"), "w+") as current:
        current.write(bedgraph_minus)

if __name__ == "__main__":
    main(sys.argv[0])
