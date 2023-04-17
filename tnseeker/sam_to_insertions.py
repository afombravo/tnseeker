import csv
import numpy as np
import os, glob
from regex import findall
from matplotlib import pyplot as plt
import argparse

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

    pathing = path_finder(folder_path)
    extractor(name_folder, folder_path, pathing, paired_ended,barcode,read_threshold,read_cut,map_quality_threshold)

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
    
    def __init__(self, contig=None, local=None, border=None, orientation=None,count=None,barcode=None,mapQ=0):
        
        self.contig = contig
        self.local = local
        self.orientation = orientation
        self.seq = border
        self.count = count
        self.mapQ = mapQ
        self.barcode = barcode or {}

def insertion_counter(current, contig, local, orientation, border):

    compiled_local = contig + local + orientation
    
    if compiled_local not in current: 
        current[compiled_local] = Insertion(contig, local, orientation, border, 1)
        
    else: 
        current[compiled_local].count += 1

    return current

def read_count(dictionary):
    count = 0

    for key in dictionary:
        count += dictionary[key].count
        
    return count

def extractor(name_folder, folder_path, pathing, paired_ended,barcode,read_threshold,read_cut,map_quality_threshold = 42):
    
    def barcode_finder():
        read = ""
        barcode = {}
        with open(f"{folder_path}/barcodes_1.txt") as current:
            for line in current:
                if "@" not in line:
                    if len(line[:-1]) != 0:
                        bar = line[:-1]
                else:
                    read = line[:-1].split(" ")[0][1:]
                    barcode[read] = bar
        return barcode

    genome_size, aligned_reads, aligned_valid_reads = 0, 0, 0
    insertion_count = {}
    
    flag_list = [0, 16]
    if paired_ended=="PE":
        flag_list = [83, 99] #[16] for single ended data #99 and 83 means that the read is the first in pair (only paired ended reads are considered as valid)

    if barcode:
        redundancy_barcode = {}
        read_to_barcode = barcode_finder()
    
    for file in pathing:
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

                            local=str(int(local)+len(sequence)-clipped) #t

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
                            if sam[0] in read_to_barcode:
                                bar = read_to_barcode[sam[0]]
                                if bar in insertion_count[key].barcode:
                                    insertion_count[key].barcode[bar] += 1
                                else:
                                    insertion_count[key].barcode[bar] = 1
                                    
                            key = (contig, local, orientation, bar)
                            if key in redundancy_barcode:
                                redundancy_barcode[key] = False
                            else:
                                redundancy_barcode[key] = True
                                    
                elif sam[0][0] == "@": #returns the genome size
                    sam = line.split('\t')
                    if "LN" in sam[-1]:
                        bp_pos = sam[-1].find(":") + 1
                        bp = int(sam[-1][bp_pos:-1])
                        genome_size = genome_size + bp

    if read_threshold:
        insertion_count,redundancy_barcode=dict_filter(insertion_count,read_cut,redundancy_barcode,barcode)
    if barcode:
        insert_parser(insertion_count,barcode,name_folder, folder_path,redundancy_barcode)
    dictionary_parser(insertion_count,barcode,folder_path,name_folder)
        
    q = plotter(insertion_count, f"Unique insertions_{name_folder}", folder_path)

    reads = f"Total Aligned Reads: {aligned_reads}\nTotal Quality Passed Reads: {aligned_valid_reads}\nFiltered Vs. Raw Read % ratio: {round(aligned_valid_reads/aligned_reads*100,2)}%\n"
    
    e = "Number of total unique insertions: {}\n".format(len(insertion_count))
    count = read_count(insertion_count)
    f = f"Total number of reads attributed to unique insertions: {count} ({round(count/aligned_valid_reads*100,2)}%)\n"
    print(reads,e,f)
    
    if barcode:
        os.remove(f"{folder_path}/barcodes_1.txt")
    
    with open("{}/library_stats_{}.txt".format(folder_path,name_folder), "w+") as current:
        current.write(reads+q+e+f)

def dict_filter(dictionary,read_cut,redundancy_barcode,barcode):
    for key in list(dictionary):
        if dictionary[key].count < read_cut:
            del dictionary[key]
            if barcode:
                for key1 in list(redundancy_barcode):
                    if (key[0] == key1[0]) & (key[1] == key1[1]) & (key[2] == key1[2]):
                        del redundancy_barcode[key1]
    return dictionary,redundancy_barcode

def insert_parser(insertion_count,barcode,name_folder, folder_path,redundancy_barcode):    
    insertions = []
    for key in insertion_count: #average bowtie Q
        insertion_count[key].mapQ = insertion_count[key].mapQ / insertion_count[key].count
        
        contig = [insertion_count[key].contig]
        local = [insertion_count[key].local]
        orientation = [insertion_count[key].orientation]
        count = [insertion_count[key].count]
        mapq = [insertion_count[key].mapQ]
        
        barcodes,reads = '',0
        if barcode:
            for bar,read in insertion_count[key].barcode.items():
                unique = 0
                key1 = (contig[0], local[0], orientation[0], bar)
                if redundancy_barcode[key1]:
                    unique = 1
                barcodes += f'{bar}:{read}:{unique};'
                reads += read
            
        insertions.append(contig + local + orientation + count + mapq + [len(insertion_count[key].barcode)] + [reads] + [barcodes])
    
    insertions.insert(0, ["#Contig"] + ["position"] + ["Orientation"] + ["Total Reads"] + \
                      ["Average MapQ"] + ["Number of different barcodes in coordinate"] + \
                      ["Total barcode Reads"] + ["Barcodes (barcode:read:unique)"])
    
    name = "barcoded_insertions_%s" % name_folder
    output_file_path = os.path.join(folder_path, name + ".csv")
    with open(output_file_path, "w", newline='') as output:
        writer = csv.writer(output)
        writer.writerows(insertions)
                
def dictionary_parser(dictionary,barcode,folder_path,name_folder):
    
    insertions = []
    
    for key in dictionary:
        
        contig = [dictionary[key].contig]
        local = [dictionary[key].local]
        orientation = [dictionary[key].orientation]
        border = [dictionary[key].seq]
        count = [dictionary[key].count]

        insertions.append(contig + local + orientation + border + count)
    insertions.insert(0, ["#Contig"] + ["position"] + ["Orientation"] + ["Transposon Border Sequence"] + ["Read Counts"])
    
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

    args = parser.parse_args()
    
    main([args.folder_path, args.name_folder, args.paired_ended, args.read_threshold, args.read_cut, args.barcode, args.map_quality_threshold])
