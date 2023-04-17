import os
import multiprocessing
import sys
from numba import njit
import numpy as np
import gzip

""" This script is for processing and trimming high-throughput sequencing data. 
    It takes as input a fastq file, a folder path to store the output, the 
    sequence of the transposon used in the experiment, and some optional 
    parameters such as whether to consider barcode information and the Phred 
    score quality threshold.

    The script has several functions that perform different tasks:
    
    1. seq2bin: converts a string to binary, and then to a numpy array in int8 format
    2. binary_subtract: used for matching two sequences based on the allowed mismatches
    3. imperfect_find: matches two sequences based on the allowed mismatches, 
    used for sequencing searching a start/end place in a read
    4. write: writes the list of processed reads to a file
    5. barcodeID: extracts the barcode from a read
    6. read_trimer: trims reads based on the presence of transposon sequences 
    and quality score
    7. cpu: returns the number of available CPU cores
    8. extractor: the main function for trimming reads and writing the output
    9. paired_ended_rearrange: rearranges paired-end reads into proper pairs
    10. main: the main function that takes the inputs, calls the functions, and writes the output
    
    The script processes the fastq file in parallel by dividing it into chunks 
    and using multiple CPU cores to trim the reads. The processed reads and 
    barcode information, if specified, are then written to files in the specified folder. 
    The script also has the capability to process paired-end reads and rearrange 
    them into proper pairs."""

def seq2bin(sequence):
    
    """ Converts a string to binary, and then to 
    a numpy array in int8 format"""

    return np.array(bytearray(sequence, 'utf8'), dtype=np.int8)

@njit
def binary_subtract(array1,array2,mismatch):
    
    """ Used for matching 2 sequences based on the allowed mismatches.
    Requires the sequences to be in numerical form"""
    
    miss=0
    for arr1,arr2 in zip(array1,array2):
        if arr1-arr2 != 0:
            miss += 1
        if miss>mismatch:
            return 0
    return 1

@njit
def imperfect_find(read,seq,mismatch): 
    
    """ Matches 2 sequences (after converting to int8 format)
    based on the allowed mismatches. Used for sequencing searching
    a start/end place in a read"""
    
    s=seq.size
    r=read.size
    fall_over_index = r-s-1
    for i,bp in enumerate(read): #range doesnt exist in njit
        comparison = read[i:s+i]
        finder = binary_subtract(seq,comparison,mismatch)
        if i > fall_over_index:
            return
        if finder != 0:
            return i

def write(listing, name, folder_path):
    text_out = folder_path + name
    if not os.path.isfile(text_out):
        with open(text_out, "w") as text_file:
            text_file.close()
    with open(text_out, "a+") as text_file:
        for item in listing:
            for item1 in item:
                text_file.write(item1 + "\n")

def barcodeID(sequence,sequence_bin,borders):
    border_up = imperfect_find(sequence_bin,borders[0],1)
    border_down = imperfect_find(sequence_bin,borders[1],1)
    if (border_up is not None) & (border_down is not None):
        return sequence[border_up+len(borders[0]):border_down]

def read_trimer(reading,sequences,quality_set,borders="",barcode_allow=False):
    processed_read = []
    barcode_pool = []
    for read in reading:
        sequence = str(read[1],"utf-8")
        sequence_bin = seq2bin(sequence)
        quality = str(read[3],"utf-8")

        border_find = imperfect_find(sequence_bin,sequences,2) #2 mismatches in Tn search sequence
        if border_find is not None:
            read[1] = sequence[border_find+len(sequences):]
            read[3] = quality[border_find+len(sequences):]

            if (len(quality_set.intersection(quality)) == 0):
                read[0],read[2] = str(read[0],"utf-8"),str(read[2],"utf-8")
                processed_read.append(read)
                if barcode_allow:
                    barcode = barcodeID(sequence,sequence_bin,borders)
                    if barcode is not None:
                        barcode_pool.append([barcode]+[read[0]])

    return [processed_read,barcode_pool]
                
def cpu():
    c = multiprocessing.cpu_count()
    if c >= 2:
        c -= 1
    pool = multiprocessing.Pool(processes = c)
    return pool, c
      
def extractor(fastq,folder_path,sequences,barcode,phred = 1):
    transposon_seq = seq2bin(sequences)
    borders = [seq2bin("CGAGGTCTCT"),seq2bin("CGTACGCT")]
    reading = []
    read_bucket=[]
    pool, cpus = cpu()
    divider = 25000
    count_total=0
    count_trimed=0
    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI" #Phred score
    if phred < 1:
        phred = 1
    quality_set = set(quality_list[:phred-1])
    with gzip.open(fastq, "rb") as current:
        for line in current:
            reading.append(line[:-1])
            
            if len(reading) == 4: #a read always has 4 lines
                read_bucket.append(reading)
                reading=[]
                count_total+=1
                
            if len(read_bucket)>=cpus*divider:
                result_objs,subdivied = [],[]
                pool, cpus = cpu()
                z, spliter_mid=0,divider
                for i in range(cpus):
                    subdivied.append(read_bucket[z:spliter_mid])
                    z += divider
                    spliter_mid += divider
                    
                for read in subdivied:
                    if not barcode:
                        result=pool.apply_async(read_trimer, args=((read,transposon_seq,quality_set)))
                    else:
                        result=pool.apply_async(read_trimer, args=((read,transposon_seq,quality_set,borders,True)))
                    result_objs.append(result)
                pool.close()
                pool.join()
                    
                result = [result.get() for result in result_objs]
                if not barcode:
                    for trimmed,barcodes in result:
                        count_trimed+=len(trimmed)
                        write(trimmed, "/processed_reads_1.fastq", folder_path)
                        
                else:
                    for trimmed,barcodes in result:
                        count_trimed+=len(trimmed)
                        write(trimmed, "/processed_reads_1.fastq", folder_path)
                        write(barcodes, "/barcodes_1.txt", folder_path)
                read_bucket = []
                
    if not barcode:
        trimmed,barcodes=read_trimer(read_bucket,transposon_seq,quality_set)
        write(trimmed, "/processed_reads_1.fastq", folder_path)
    else:
        trimmed,barcodes=read_trimer(read_bucket,transposon_seq,quality_set,borders)
        write(trimmed, "/processed_reads_1.fastq", folder_path)
        write(barcodes, "/barcodes_1.txt", folder_path)
    count_trimed+=len(trimmed)
    
    text_out = folder_path + "/trimming_log.log"
    with open(text_out, "w+") as text_file:
        text_file.write(f"Total reads trimmed: {count_trimed}\nTotal reads in file: {count_total}\nPercent of passing reads: {count_trimed/count_total*100}\n")
    read_bucket,result = [],[]

def paired_ended_rearrange(fastq2,folder_path):
    names=set()
    duplicated=set()
    with open(folder_path+"/processed_reads_1.fastq") as current:
        for line in current:
            if "@" in line:
                title = line.split(" ")[0]
                title = title.split(":")
                title = title[3]+title[4]+title[5]+title[6]
                if title not in names:
                    names.add(title)
                else:
                    duplicated.add(title)
    
    reading,read_bucket=[],[]
    with gzip.open(fastq2, "rb") as current:
        for line in current:
            reading.append(line[:-1])
            if len(reading) == 4: #a read always has 4 lines
                reading[0]=str(reading[0],"utf-8")
                title = reading[0].split(" ")[0]
                title = title.split(":")
                title = title[3]+title[4]+title[5]+title[6]
                if ("@" in reading[0]) & ((title in names) or (title in duplicated)):
                    if title in duplicated:
                        duplicated.remove(title)
                    else:
                        names.remove(title)
                    reading[1]=str(reading[1],"utf-8")
                    reading[2]=str(reading[2],"utf-8")
                    reading[3]=str(reading[3],"utf-8")
                    read_bucket.append(reading)
                    if len(read_bucket)>1000000:
                        write(read_bucket, "/processed_reads_2.fastq", folder_path,False)
                        read_bucket=[]
                reading=[]
                        
    write(read_bucket, "/processed_reads_2.fastq", folder_path,False)
    
def main(argv):

    fastq1=argv[0] #"C:/Users/afons/Documents/Tn-Seq/HiSeq5_data/NT12004/NT12004_HiSeq5/HFVHGBBXY_Milion_HiSeq_19s005084-1-1_Bravo_lane6sample2_sequence.fastq"
    folder_path = argv[1] #"C:/Users/afons/Documents/Tn-Seq/Carlos/liquid"
    sequences = argv[2]
    paired = argv[3]
    phred = int(argv[5])
    
    barcode = False
    if argv[4] == "True":
        barcode = True

    print("Trimming Sequences")
    try:
        extractor(fastq1,folder_path,sequences,barcode,phred)
    except Exception as e:
        print(e)
    
    if paired == "PE":
        fastq2=argv[6]
        paired_ended_rearrange(fastq2,folder_path)
            
if __name__ == "__main__":
    if len(sys.argv) > 6:
        argv = sys.argv[1:7] if sys.argv[4] == "PE" else sys.argv[1:6]
    main(argv)
