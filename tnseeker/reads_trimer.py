import os
import multiprocessing
import sys
import gzip
from fast2q.fast2q import border_finder,seq2bin
from tnseeker.extras.helper_functions import colourful_errors,file_finder,variables_parser

""" This script is for processing and trimming fastq files. 
    It takes as input a fastq file, a folder path to store the output, 
    and some optional parameters such as whether to consider barcode 
    information and the Phred score quality threshold."""

def write(listing, name, folder_path):
    text_out = folder_path + name
    if not os.path.isfile(text_out):
        with open(text_out, "w") as text_file:
            text_file.close()
    with open(text_out, "a+") as text_file:
        for item in listing:
            for item1 in item:
                text_file.write(item1 + "\n")

def barcode_finder(sequence,sequence_bin,variables):
    border_up = border_finder(variables["borders_bin"][0],sequence_bin,variables['barcode_up_miss'])
    if border_up is not None:
        start_place = border_up+len(variables["borders_bin"][0])
        border_down = border_finder(variables["borders_bin"][1],sequence_bin,variables['barcode_down_miss'],start_place)
        if border_down is not None:
            return sequence[start_place:border_down]
    return None

def read_trimer(read_chunk):
    
    processed_read = []
    for read in read_chunk:
        sequence = str(read[1],"utf-8")
        sequence_bin = seq2bin(sequence)
        quality = str(read[3],"utf-8")

        border_find = border_finder(variables["sequence_bin"],sequence_bin,variables["tn_mismatches"]) 
        if border_find is not None:
            
            start = border_find+len(variables["sequence_bin"])
            if variables['trimmed_after_tn'] != -1:
                end = start+variables['trimmed_after_tn']
                read[1] = sequence[start:end]
                read[3] = quality[start:end]
            else:
                read[1] = sequence[start:]
                read[3] = quality[start:]

            if (len(variables["quality_set"].intersection(quality)) == 0):
                read[0],read[2] = str(read[0],"utf-8"),str(read[2],"utf-8")
                read[0] = read[0].split(" ")[0]
                if variables["barcode"]:
                    barcode = barcode_finder(sequence,sequence_bin,variables)
                    if barcode is not None:
                        if (len(variables["quality_set_bar_up"].intersection(quality)) == 0) & (len(variables["quality_set_bar_down"].intersection(quality)) == 0):
                            read[0] += f":BC:{barcode}"
                processed_read.append(read)

    return processed_read
                
def extractor():
    
    variables["sequence_bin"] = seq2bin(variables["sequence"])  
    reading = []
    read_bucket=[]
    divider = 250000
    count_total=0
    count_trimed=0
    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI" #Phred score
    if variables["phred"] < 1:
        variables["phred"] = 1
    variables["quality_set"] = set(quality_list[:variables["phred"]-1])
    
    if variables["barcode"]:
        variables["borders_bin"] = [seq2bin(variables["barcode_up"]),seq2bin(variables["barcode_down"])]
        variables["quality_set_bar_up"] = set(quality_list[:variables["barcode_up_phred"]-1])
        variables["quality_set_bar_down"] = set(quality_list[:variables["barcode_down_phred"]-1])

    file_number = len(variables['fastq_file'])
    file_counter = 0
    for file in variables['fastq_file']:
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
                        
                    if len(read_bucket)>=variables['cpus']*divider:
                        pool = multiprocessing.Pool(processes = variables['cpus'])
                        result_objs,subdivied = [],[]
                        z, spliter_mid=0,divider
                        for i in range(variables['cpus']):
                            subdivied.append(read_bucket[z:spliter_mid])
                            z += divider
                            spliter_mid += divider

                        for read_chunk in subdivied:
                            result=pool.apply_async(read_trimer,args=(read_chunk,))
                            result_objs.append(result)
                        pool.close()
                        pool.join()
                            
                        result = [result.get() for result in result_objs]
                        for trimmed in result:
                            count_trimed+=len(trimmed)
                            write(trimmed, "/processed_reads_1.fastq", variables['directory'])
                        read_bucket = []
                        
        except Exception as e:
            colourful_errors("WARNING",f'Error parsing {file} due to {e}')
            
    trimmed=read_trimer(read_bucket)
    write(trimmed, "/processed_reads_1.fastq", variables['directory'])

    count_trimed+=len(trimmed)
    text_out = variables['directory'] + "/trimming_log.log"
    with open(text_out, "w+") as text_file:
        text_file.write(f"""    ####\nTN TRIMMING INFO\n    ####\n\nTotal reads trimmed: {count_trimed}\nTotal reads in file: {count_total}\nPassing reads: {round(count_trimed/count_total*100,2)}%\n\n    ####\nBOWTIE2 INFO\n    ####\n\n""")

def paired_ended_rearrange():
    fastq2=file_finder(variables['sequencing_files_r'],['*.gz'])
    names=set()
    duplicated=set()
    reading = []
    with open(variables['directory']+"/processed_reads_1.fastq") as current:
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
                            write(read_bucket, "/processed_reads_2.fastq", variables['directory'])
                            read_bucket=[]
                    reading=[]
                        
    write(read_bucket, "/processed_reads_2.fastq", variables['directory'])

def main(argv):

    global variables
    variables = variables_parser(argv)
    variables['fastq_file']=file_finder(variables['sequencing_files'],['*.gz'])
    variables['phred'] = int(variables['phred'])
    variables['tn_mismatches'] = int(variables['tn_mismatches'])
    variables['trimmed_after_tn'] = int(variables['trimmed_after_tn'])
    variables["barcode"] = variables["barcode"] == "True"
    variables["barcode_up"] = None if variables["barcode_up"] == "None" else variables["barcode_up"]
    variables["barcode_down"] = None if variables["barcode_down"] == "None" else variables["barcode_down"]
    variables["barcode_up_miss"] = int(variables["barcode_up_miss"])
    variables["barcode_down_miss"] = int(variables["barcode_down_miss"]) 
    variables["barcode_up_phred"] = int(variables["barcode_up_phred"]) 
    variables["barcode_down_phred"] = int(variables["barcode_down_phred"]) 
    variables['cpus'] = int(variables['cpus'])
    variables['pool'] = multiprocessing.Pool(processes = variables['cpus'])

    extractor()

    if variables['seq_type'] == "PE":
        if 'sequencing_files_r' in variables:
            paired_ended_rearrange()
            
if __name__ == "__main__":
    main(sys.argv[0])
    multiprocessing.set_start_method("spawn")

