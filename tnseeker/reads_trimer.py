import os
import multiprocessing
from tnseeker.extras.helper_functions import colourful_errors
import sys
from numba import njit
import numpy as np
import gzip
import glob

""" This script is for processing and trimming high-throughput sequencing data. 
    It takes as input a fastq file, a folder path to store the output, the 
    sequence of the transposon used in the experiment, and some optional 
    parameters such as whether to consider barcode information and the Phred 
    score quality threshold."""

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
def imperfect_find(read,seq,mismatch,start_place=0): 
    
    """ Matches 2 sequences (after converting to int8 format)
    based on the allowed mismatches. Used for sequencing searching
    a start/end place in a read"""
    
    s=seq.size
    r=read.size
    fall_over_index = r-s-1
    for i,bp in enumerate(read[start_place:]): 
        comparison = read[start_place+i:s+start_place+i]
        finder = binary_subtract(seq,comparison,mismatch)
        if i > fall_over_index:
            return None
        if finder != 0:
            return i+start_place

def write(listing, name, folder_path):
    text_out = folder_path + name
    if not os.path.isfile(text_out):
        with open(text_out, "w") as text_file:
            text_file.close()
    with open(text_out, "a+") as text_file:
        for item in listing:
            for item1 in item:
                text_file.write(item1 + "\n")

def barcodeID(sequence,sequence_bin,borders,miss_up,miss_down):
    border_up = imperfect_find(sequence_bin,borders[0],miss_up)
    if border_up is not None:
        start_place = border_up+len(borders[0])
        border_down = imperfect_find(sequence_bin,borders[1],miss_down,start_place)
        if border_down is not None:
            return sequence[start_place:border_down]
    return None

def read_trimer(reading,sequences,quality_set,mismatches,trimming_len,miss_up,\
                miss_down,quality_set_bar_up,quality_set_bar_down,borders="",barcode_allow=False):
    
    processed_read = []
    for read in reading:
        sequence = str(read[1],"utf-8")
        sequence_bin = seq2bin(sequence)
        quality = str(read[3],"utf-8")

        border_find = imperfect_find(sequence_bin,sequences,mismatches) 
        if border_find is not None:
            
            start = border_find+len(sequences)
            if trimming_len != -1:
                end = start+trimming_len
                read[1] = sequence[start:end]
                read[3] = quality[start:end]
            else:
                read[1] = sequence[start:]
                read[3] = quality[start:]

            if (len(quality_set.intersection(quality)) == 0):
                read[0],read[2] = str(read[0],"utf-8"),str(read[2],"utf-8")
                read[0] = read[0].split(" ")[0]
                if barcode_allow:
                    barcode = barcodeID(sequence,sequence_bin,borders,miss_up,miss_down)
                    if barcode is not None:
                        if (len(quality_set_bar_up.intersection(quality)) == 0) & (len(quality_set_bar_down.intersection(quality)) == 0):
                            read[0] += f":BC:{barcode}"
                processed_read.append(read)

    return processed_read
                
def extractor(fastq,folder_path,sequences,barcode,barcode_upstream,barcode_downstream,\
              mismatches,trimming_len,miss_up,miss_down,phred_up,phred_down,
              cpus,pool,phred = 1):
    
    transposon_seq = seq2bin(sequences)  
    reading = []
    read_bucket=[]
    divider = 250000
    count_total=0
    count_trimed=0
    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI" #Phred score
    if phred < 1:
        phred = 1
    quality_set = set(quality_list[:phred-1])
    
    quality_set_bar_up,quality_set_bar_down = None,None
    if barcode:
        borders = [seq2bin(barcode_upstream),seq2bin(barcode_downstream)]
        quality_set_bar_up = set(quality_list[:phred_up-1])
        quality_set_bar_down = set(quality_list[:phred_down-1])

    file_number = len(fastq)
    file_counter = 0
    for file in fastq:
        file_counter += 1
        colourful_errors("INFO",
            f"Processing {file_counter} out of {file_number} fastq files.")
        try:
            with gzip.open(file, "rb") as current:
                for line in current:
                    reading.append(line[:-1])
                    
                    if len(reading) == 4: 
                        read_bucket.append(reading)
                        reading=[]
                        count_total+=1
                        
                    if len(read_bucket)>=cpus*divider:
                        pool = multiprocessing.Pool(processes = cpus)
                        result_objs,subdivied = [],[]
                        z, spliter_mid=0,divider
                        for i in range(cpus):
                            subdivied.append(read_bucket[z:spliter_mid])
                            z += divider
                            spliter_mid += divider

                        for read in subdivied:
                            if not barcode:
                                result=pool.apply_async(read_trimer, 
                                                        args=((read,transposon_seq,quality_set,mismatches,trimming_len,miss_up,miss_down,\
                                                               quality_set_bar_up,quality_set_bar_down)))
                            else:
                                result=pool.apply_async(read_trimer, 
                                                        args=((read,transposon_seq,quality_set,mismatches,trimming_len,miss_up,miss_down,\
                                                               quality_set_bar_up,quality_set_bar_down,borders,True)))
                            result_objs.append(result)
                        pool.close()
                        pool.join()
                            
                        result = [result.get() for result in result_objs]
                        for trimmed in result:
                            count_trimed+=len(trimmed)
                            write(trimmed, "/processed_reads_1.fastq", folder_path)
 
                        read_bucket = []
                        
        except Exception:
            colourful_errors("WARNING",
                f'Error parsing {file}')
            

    trimmed=read_trimer(read_bucket,transposon_seq,quality_set,mismatches,trimming_len,miss_up,miss_down,
                                 quality_set_bar_up,quality_set_bar_down)
    write(trimmed, "/processed_reads_1.fastq", folder_path)

    count_trimed+=len(trimmed)
    
    text_out = folder_path + "/trimming_log.log"
    with open(text_out, "w+") as text_file:
        text_file.write(f"Total reads trimmed: {count_trimed}\nTotal reads in file: {count_total}\nPercent of passing reads: {count_trimed/count_total*100}\n")
    read_bucket,result = [],[]

def paired_ended_rearrange(fastq2,folder_path):
    names=set()
    duplicated=set()
    reading = []
    with open(folder_path+"/processed_reads_1.fastq") as current:
        for line in current:
            reading.append(line[:-1])
            if len(reading) == 4: 
                title = reading[0].split(" ")[0]
                title = title.split(":")
                title = title[3]+title[4]+title[5]+title[6]
                if title not in names:
                    names.add(title)
                else:
                    duplicated.add(title)
                reading=[]
    
    reading,read_bucket=[],[]
    for file in fastq2:
        with gzip.open(file, "rb") as current:
            for line in current:
                reading.append(line[:-1])
                if len(reading) == 4: 
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
                            write(read_bucket, "/processed_reads_2.fastq", folder_path)
                            read_bucket=[]
                    reading=[]
                        
    write(read_bucket, "/processed_reads_2.fastq", folder_path)
   
def folder_sequence_parser(folder):
    pathing = []
    for exten in ['*.gz']:
        for filename in glob.glob(os.path.join(folder, exten)):
            pathing.append(filename) 
    return pathing

def main(argv):

    fastq1=folder_sequence_parser(argv[0])
    folder_path = argv[1] 
    sequences = argv[2]
    paired = argv[3]
    phred = int(argv[5])
    mismatches = int(argv[-3])
    trimming_len = int(argv[-2])
    cpus = int(argv[-1])
    pool = multiprocessing.Pool(processes = cpus)
    
    barcode,barcode_upstream,barcode_downstream,miss_up,miss_down,phred_up,phred_down = False,None,None,None,None,None,None
    if argv[4] == "True":
        barcode = True
        barcode_upstream = argv[-9]
        barcode_downstream = argv[-8]
        miss_up = int(argv[-7])
        miss_down = int(argv[-6])
        phred_up = int(argv[-5])
        phred_down = int(argv[-4])

    try:
        extractor(fastq1,folder_path,sequences,barcode,barcode_upstream,\
                  barcode_downstream,mismatches,trimming_len,miss_up,miss_down,\
                  phred_up,phred_down,cpus,pool,phred)
    except Exception as e:
        print(e)
    
    if paired == "PE":
        fastq2=folder_sequence_parser(argv[6])
        paired_ended_rearrange(fastq2,folder_path)
            
if __name__ == "__main__":
    if len(sys.argv) > 7:
        argv = sys.argv[1:13] if sys.argv[4] == "PE" else sys.argv[1:12]
        if argv[4] == "True":
            argv.append(sys.argv[-7],sys.argv[-6],sys.argv[-5])
    main(argv)
    
    multiprocessing.set_start_method("spawn")

