import statsmodels.stats.multitest
from Bio import SeqIO
import os
import glob
import numpy as np
import scipy
import re
from tnseeker.extras.possion_binom import PoiBin
from scipy.stats import binomtest
import multiprocessing
import matplotlib.pyplot as plt
import sys
from numba import njit
import pkg_resources
import pandas as pd
import csv
import subprocess
from pathlib import Path
from colorama import Fore
import datetime
import matplotlib
matplotlib.use('Agg')

"""
Disclaimer: 
    I am aware this script is messy and complex. 
    It was 5y in the making and changed a lot over time to accommodate endless 
    features, exceptions, and use cases.
    It´s more of a pathwork than anything else at the moment.
    Please don´t judge me too harshly.
"""

def inputs(argv):
    ''' The inputs function initializes a global Variables instance with various attributes 
    from a given list of command-line arguments. It sets up parameters for directory, 
    strain, annotation, subdomain length, p-value, and output, as well as true positives 
    and true negatives file paths. '''

    global variables
    variables = Variables()
    variables.directory = argv[0]  # folder with the unique insertions file
    variables.strain = argv[1]  # strain name, and annotation file name
    variables.annotation_type = argv[2]
    variables.annotation_folder = argv[3]
    variables.subdomain_length = [float(argv[4]), float(argv[5])] #upstream, downstream
    variables.pvalue = float(argv[6])
    variables.ir_size_cutoff = int(argv[7])
    variables.output_name = variables.strain + "_alldomains"
    variables.domain_uncertain_threshold = float(argv[8])  # 0.75
    variables.biggest_gene = 0

    variables.true_positives = pkg_resources.resource_filename(
        __name__, 'data/true_positives.fasta')
    variables.true_negatives = pkg_resources.resource_filename(
        __name__, 'data/true_negatives.fasta')

def path_finder():
    ''' The path_finder function searches for and sets file paths for both insertions file and the 
    correct annotation (gb/gff) based on the global variables instance. It uses the nested sub_path_finder 
    function to search for specific file types and names within specified folders. '''

    def sub_path_finder(folder, extenction, search):
        for filename in glob.glob(os.path.join(folder, extenction)):
            test1 = filename.find(search)
            if test1 != -1:
                return filename

    variables.insertion_file_path = sub_path_finder(
        variables.directory, '*.csv', "all_insertions")

    if variables.annotation_type == "gb":
        extention = '*.gb'
        if sub_path_finder(variables.annotation_folder, extention, variables.strain) == None:
            extention = '*.gbk'
        variables.annotation_file_paths = [sub_path_finder(
            variables.annotation_folder, extention, variables.strain)]

    elif variables.annotation_type == "gff":
        variables.annotation_file_paths = [sub_path_finder(
            variables.annotation_folder, '*.gff', variables.strain)]
        variables.annotation_file_paths.append(sub_path_finder(
            variables.annotation_folder, '*.fasta', variables.strain))

        if len(variables.annotation_file_paths) < 2:
            print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.RED}FATAL{Fore.RESET}] Check there is a .fasta and .gff file in the indicated folder")
            raise Exception


def output_writer(output_folder, name_folder, output_file):
    ''' The output_writer function writes the given output_file data as a CSV file 
    to the specified output_folder with the provided name_folder. 
    It ensures proper formatting and line handling by using the csv module. '''

    output_file_path = os.path.join(output_folder, name_folder + ".csv")
    with open(output_file_path, "w", newline='') as output:
        writer = csv.writer(output)
        writer.writerows(output_file)


class Variables():

    ''' The Variables class is a comprehensive container for various attributes 
    and methods related to genomic data analysis. It initializes a range of attributes, 
    such as directory paths, strain information, annotation details, contig information, 
    genome sequence and length, and statistical parameters, among others. 
    Users can provide custom values for these attributes, or the class will 
    assign default values when initialized. 

    The class also includes two methods for handling motif-related calculations:

        motif_compiler: This method generates a list of dinucleotide motifs (di_motivs) 
        and their corresponding compiled regular expressions (regex_compile) for 
        pattern matching. It covers all possible combinations of the four 
        DNA bases (A, T, C, G) in dinucleotide motifs.

        normalizer: This method calculates the probability of having a transposon 
        insertion in each dinucleotide motif by counting the entire motif content 
        of the genome and dividing the transposon motif count by the motif 
        genome count. It stores the results in the chance_motif_tn attribute. 
        If any entry in the resulting array is greater than or equal to 1, 
        it is set to 0.999999. This is related with downstream calculations using
        the poibin module.

    Overall, the Variables class serves as a central storage and management unit 
    for various attributes and calculations essential for genomic data analysis. 
    It streamlines data handling, allowing for smoother execution of analysis tasks. '''

    def motif_compiler(self):
        dna = ["A", "T", "C", "G"]
        prog, di_motivs = [], []
        for letter1 in dna:
            for letter2 in dna:
                di_motivs.append(f"{letter1}{letter2}")
                motiv = f"(?=({letter1}{letter2}))"
                prog.append(re.compile(motiv, re.IGNORECASE))
        return prog, di_motivs

    def __init__(self, directory=None, strain=None, annotation_type=None, annotation_folder=None, pan_annotation=None,
                 output_name=None, true_positives=None, true_negatives=None,
                 insertion_file_path=None, annotation_file_paths=None, borders_contig={}, orientation_contig={},
                 insertions_contig={}, genome_seq={}, genome_length=0, annotation_contig=None, total_insertions=None,
                 positive_strand_tn_ratio=None, transposon_motiv_count=None, transposon_motiv_freq=None,
                 chance_motif_tn=None, orientation_contig_plus=None, orientation_contig_neg=None,
                 subdomain_length=None, pvalue=None):

        self.directory = directory
        self.strain = strain
        self.annotation_type = annotation_type
        self.annotation_folder = annotation_folder
        self.pan_annotation = pan_annotation
        self.output_name = output_name
        self.true_positives = true_positives
        self.true_negatives = true_negatives
        self.domain_iteration = [1, 1.5, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                                 1024, 2048, 4096, 8192, 16384, 32768]  # """ make user defined """
        self.subdomain_length = subdomain_length or [0.1, 0.9]
        self.pvalue = pvalue or 0.1
        self.insertion_file_path = insertion_file_path
        self.annotation_file_paths = annotation_file_paths
        self.borders_contig = borders_contig
        self.orientation_contig = orientation_contig
        self.insertions_contig = insertions_contig
        self.genome_seq = genome_seq
        self.genome_length = genome_length
        self.annotation_contig = annotation_contig
        self.total_insertions = total_insertions
        self.positive_strand_tn_ratio = positive_strand_tn_ratio
        self.transposon_motiv_count = transposon_motiv_count or dict()
        self.transposon_motiv_freq = transposon_motiv_freq or dict()
        self.chance_motif_tn = chance_motif_tn or dict()
        self.orientation_contig_plus = orientation_contig_plus or dict()
        self.orientation_contig_neg = orientation_contig_neg or dict()
        self.regex_compile, self.di_motivs = self.motif_compiler()

    def normalizer(self,contig):
        motif_genome = np.zeros(16)
        # Counts the entire motifv content of the genome; determining the probability of having an insertion in each motif
        atcg = count_GC([self.genome_seq[contig]], self.regex_compile)
        for i, element in enumerate(motif_genome):
            motif_genome[i] = element + atcg[i]

        # probability oh having a motif with a transposon in it
        self.chance_motif_tn[contig] = np.divide(
            self.transposon_motiv_count[contig], motif_genome)
        for i, entry in enumerate(self.chance_motif_tn[contig]):
            if entry >= 1:
                self.chance_motif_tn[contig][i] = 0.999999
        return np.array([[n] for n in self.chance_motif_tn[contig]])


class Significant:

    ''' The Significant class is a container for storing information about a 
    gene or domain. It initializes attributes such as domain insertions, 
    p-value, domain length, orientation ratio, orientation p-value, 
    gene region, and essentiality. Used for the final processing of essentiality'''

    def __init__(self, dom_insert=None, pvalue=None, dom_len=None, ratio_orient=None,
                 orient_pvalue=None, domain_part=None, essentiality="Non-Essential"):
        self.dom_insert = dom_insert
        self.pvalue = pvalue
        self.dom_len = dom_len
        self.ratio_orient = ratio_orient
        self.orient_pvalue = orient_pvalue
        self.domain_part = domain_part
        self.essentiality = essentiality


class Gene:

    ''' The Gene class is a container for storing information about a gene, 
    including its attributes such as start, end, orientation, domains, identity, 
    product, contig, and various data related to insertions, GC content, 
    and domain significance. It initializes these attributes with provided or 
    default values. Used as the first step for storing all processed gene information 
    from the annotation files'''

    def __init__(self, gene=None, start=None, end=None, orientation=None, domains=None, identity=None,
                 product=None, contig=None, matrix=None, domain_insertions_total=None, GC_content=None,
                 domain_notes=None, motif_seq=None, subdomain_insert_orient_plus=None,
                 subdomain_insert_orient_neg=None, significant=None):

        self.gene = gene
        self.start = start
        self.end = end
        self.length = end - start
        self.orientation = orientation
        self.domains = domains
        self.identity = identity
        self.product = product
        self.contig = contig
        self.gene_insert_matrix = matrix
        self.domain_insertions_total = domain_insertions_total
        self.GC_content = GC_content
        self.domain_notes = domain_notes
        self.motif_seq = motif_seq
        self.subdomain_insert_orient_plus = subdomain_insert_orient_plus or dict()
        self.subdomain_insert_orient_neg = subdomain_insert_orient_neg or dict()
        self.significant = significant or dict()

def subprocess_cmd(command):
    try:
        return subprocess.check_output(command)
    except subprocess.CalledProcessError as e:
        return e.output.decode()

def blast_maker():
    
    def tblastn(file):
        
        out = f"{blast_db}/{Path(variables.strain).stem}_{Path(file).stem}.txt"

        send = ["tblastn",
                "-query",
                file,
                "-db",
                fasta,
                "-evalue",
                "0.00001",
                "-out",
                out,
                "-outfmt",
                "6"]
        
        subprocess_cmd(send)

        found_true = {}
        with open(out) as current:
            for line in current:
                line = line.split("\t")
                first = int(line[8])
                last = int(line[9])
                if first > last:
                    orientation = "-"
                    start = last
                    end = first
                else:
                    orientation = "+"
                    start = first
                    end = last
                found_true[line[0]] = Gene(gene=line[0],
                                           contig=line[1],
                                           start=start,
                                           end=end,
                                           orientation=orientation)

        return found_true
    
    blast_db = f"{variables.directory}/blast_db"
    if not os.path.isdir(blast_db):
        os.mkdir(blast_db)
    
    if variables.annotation_type == "gb":
        fasta = f"{blast_db}/{Path(variables.strain).stem}.fasta"
        with open(variables.annotation_file_paths[0], "r") as input_handle:
            with open(fasta, "w") as output_handle:
                sequences = SeqIO.parse(input_handle, "genbank")
                SeqIO.write(sequences, output_handle, "fasta")
    else:
        fasta_in = variables.annotation_file_paths[1]
        fasta = f"{blast_db}/{Path(variables.strain).stem}.fasta"
        send = ["cp",
                fasta_in,
                fasta]
        subprocess_cmd(send)

    send = ["makeblastdb",
            "-in",
            fasta,
            "-dbtype",
            "nucl"]
    
    subprocess_cmd(send)

    variables.true_positives = tblastn(variables.true_positives)
    variables.true_negatives = tblastn(variables.true_negatives)
    
    print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.GREEN}INFO{Fore.RESET}] Found {len(variables.true_positives)} benchmark essential genes.")
    
    compiled_benchmark = {**variables.true_positives, **variables.true_negatives}
    for gene in compiled_benchmark:
        lenght = compiled_benchmark[gene].end-compiled_benchmark[gene].start
        if lenght > variables.biggest_gene:
            variables.biggest_gene = lenght
    
    return compiled_benchmark

def domain_resizer(domain_size_multiplier, basket):
    ''' The domain_resizer function takes a domain size multiplier and a dictionary 
    of genes (basket) as input, and resizes gene domains by calculating the new 
    domain size based on genome length, total insertions, and domain size multiplier. 
    It then updates the domain boundaries for each gene in the basket.'''

    domain_size = int(variables.genome_length /
                      variables.total_insertions * domain_size_multiplier)
    for key in basket:
        local_stop = [basket[key].start]
        start_iterator = basket[key].start
        end = basket[key].end

        divider = int((end - start_iterator) / domain_size)

        if divider >= 1:
            for i in range(divider):
                # creates all the subdomains
                local_stop.append(start_iterator+domain_size)
                start_iterator += domain_size

        if (end - (local_stop[-1] + domain_size)) > domain_size:
            del local_stop[-1]

        if local_stop[-1] != end:
            local_stop.append(end)
        # the end of the gene is always required in duplicate
        local_stop.append(end)
        basket[key].domains = local_stop
    return basket


def gene_insertion_matrix(basket):
    ''' The gene_insertion_matrix function generates gene insertion matrices for 
    each gene in the given basket by considering genome sequence, insertions, 
    borders, and orientations of transposons. It updates the basket with the 
    calculated gene insertion matrix for each gene.'''

    def contig_matrix_generator(orient, genome_orient_plus_matrix, genome_orient_neg_matrix,
                                borders, genome_borders_matrix, inserts, genome_insert_matrix, contig):

        for j, local in enumerate(inserts):
            local = int(local-1)
            genome_insert_matrix[local] = 1
            genome_borders_matrix[local] = borders[j]
            if orient[j] == "+":
                genome_orient_plus_matrix[local] = 1
            else:
                genome_orient_neg_matrix[local] = 1

        return genome_insert_matrix, genome_borders_matrix, genome_orient_plus_matrix, genome_orient_neg_matrix

    print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.GREEN}INFO{Fore.RESET}] Compiling insertion matrix")

    for contig in variables.genome_seq:
        contig_size = len(variables.genome_seq[contig])
        genome_insert_matrix = np.zeros(contig_size, dtype=np.int8)
        genome_borders_matrix = np.zeros(contig_size, dtype='<U2')
        genome_orient_plus_matrix = np.zeros(contig_size, dtype=np.int8)
        genome_orient_neg_matrix = np.zeros(contig_size, dtype=np.int8)
        variables.orientation_contig_plus[contig] = {}
        variables.orientation_contig_neg[contig] = {}

        if len(variables.insertions_contig[contig]) != 0:
            inserts = np.array(variables.insertions_contig[contig], dtype=int)
            borders = np.array(variables.borders_contig[contig])
            orient = np.array(variables.orientation_contig[contig])
        else:
            inserts = np.zeros(contig_size, dtype=np.int8)
            borders = np.zeros(contig_size, dtype='<U2')
            orient = np.zeros(contig_size, dtype=np.int8)

        variables.insertions_contig[contig],\
            variables.borders_contig[contig],\
            variables.orientation_contig_plus[contig],\
            variables.orientation_contig_neg[contig] = \
            contig_matrix_generator(orient, genome_orient_plus_matrix, genome_orient_neg_matrix,
                                    borders, genome_borders_matrix, inserts, genome_insert_matrix, contig)

    for key in basket:
        start = basket[key].start
        end = basket[key].end
        for contig in variables.insertions_contig:
            # calculating local insertion transposon density in bp windows (size variable)
            if basket[key].contig == contig:
                basket[key].gene_insert_matrix = variables.insertions_contig[contig][start-1:end-1]

    return basket


def count_GC(seq, prog):
    ''' The count_GC function takes a sequence and a list of compiled regex patterns, 
    then counts the occurrences of each pattern in the sequence. 
    It returns a list with the counts for each pattern.'''

    result = []
    for motiv in prog:
        for insert in seq:
            result.append(len([m.start() for m in re.finditer(motiv, insert)]))
    return result


def motiv_compiler(seq, prog):
    ''' The motiv_compiler function takes a sequence and a list of compiled 
    regex patterns, then counts the occurrences of each non-empty pattern in 
    the sequence. It returns an array containing the counts for each pattern.'''

    motiv_inbox = np.zeros(16)
    for s, insertion in enumerate(seq):
        if insertion != "":
            for i, motiv in enumerate(prog):
                if [m.start() for m in re.finditer(motiv, insertion)] != []:
                    motiv_inbox[i] = motiv_inbox[i] + 1
                    break
    return motiv_inbox


def poisson_binomial(events, motiv_inbox, tn_chance):
    ''' The poisson_binomial function calculates the Poisson binomial 
    cumulative distribution function for a given set of events, motif occurrences, 
    and transposon insertion chances. It returns the absolute value of the 
    calculated p-value.'''

    def chance_matrix(events, motiv_inbox, tn_chance):
        p_effective = np.array(())
        for i, (chance, number) in enumerate(zip(tn_chance, events)):
            for i in range(number):
                p_effective = np.append(p_effective, chance)
        sucess = np.sum(motiv_inbox)
        # due to insertion redundancy (+ and - strand), sometimes there are more insertions than bp
        if sucess > np.sum(events):
            sucess = np.sum(events)
        return sucess, p_effective
    
    sucess, p_effective = chance_matrix(events, motiv_inbox, tn_chance)
    pb = PoiBin(p_effective)

    # in some cases of highly biassed transposon, essential genes can be so significant that p < 0 (probably some bug on the poisson code)
    return abs(pb.cdf(int(sucess)))


def essentials(chunk, variables):
    ''' The essentials function takes a chunk of genes and related variables as 
    input and analyzes each gene for essential domains based on various criteria 
    such as insertion patterns, orientation ratios, and statistical significance 
    (p-values). The function works as follows:

        1. Initializes data structures to store essential domain information for 
        each gene.

        2. Processes each gene in the chunk, identifying domains with no insertions 
        or specific insertion patterns.

        3. For each gene, it evaluates the ratio of insertions with positive and 
        negative orientations and calculates a p-value using a binomial test to 
        determine the significance of the ratio.

        4. If the domain spans the entire gene or if there are sub-domains within 
        the gene, it calculates essentiality metrics for each domain or sub-domain.

        5. It calls the domain_maker function to update the essentiality 
        information for each domain of the gene, such as p-value, number of insertions, 
        domain length, orientation ratio, and orientation p-value.

        6. The function iterates through all the genes in the chunk, updating 
        their essentiality information.

        7. After processing all the genes in the chunk, the function returns the 
        updated chunk with essentiality information for each gene.'''

    def ratio_insertions(p, subdomain_insert_orient_neg_cluster, subdomain_insert_orient_pos_cluster):
        total_orientation = subdomain_insert_orient_pos_cluster + \
            subdomain_insert_orient_neg_cluster
        if total_orientation != 0:
            ratio_orientation = subdomain_insert_orient_pos_cluster / \
                total_orientation  # zero means no neg insertions, only positive
            pvalue = binomtest(subdomain_insert_orient_pos_cluster,
                               total_orientation, p, alternative='two-sided').pvalue
        else:
            ratio_orientation, pvalue = "N/A", "N/A"
        return ratio_orientation, pvalue

    def domain_maker(chunk, key, domain, domain_motivs, motiv_insertions, subdomain_insert_orient_neg_cluster, subdomain_insert_orient_pos_cluster):
        ratio_orientation, orient_pvalue = ratio_insertions(variables.positive_strand_tn_ratio,
                                                            subdomain_insert_orient_neg_cluster,
                                                            subdomain_insert_orient_pos_cluster)
        chunk[key].significant[domain] = Significant()
        chunk[key].significant[domain].dom_insert = sum(motiv_insertions)
        chunk[key].significant[domain].dom_len = sum(domain_motivs)
        chunk[key].significant[domain].ratio_orient = ratio_orientation
        chunk[key].significant[domain].orient_pvalue = orient_pvalue
        chunk[key].significant[domain].domain_part = domain
        try:
            chunk[key].significant[domain].pvalue = poisson_binomial(
                domain_motivs, motiv_insertions, variables.chance_motif_tn[chunk[key].contig])  # / (1/sum(domain_motivs))
        except Exception as e:
            print(e,key)
        return chunk

    # annexes all the genes which fit the criteria of essentiality
    for key in chunk:

        chunk[key].significant = {}  # clears any overlapping dictionaries
        # total number of insertions in gene
        insertions = chunk[key].domain_insertions_total
        domains = chunk[key].domains[:-1]
        start = chunk[key].start

        if sum(insertions) == 0:  # gene has no insertions

            # number of each motif in the gene (sum is gene lenght)
            domain_motivs = np.sum(chunk[key].GC_content, axis=0)*2
            # number of each insertion per motif in the gene (sum is total insertions)
            motiv_insertions = np.sum(chunk[key].motif_seq, axis=0)

            chunk = domain_maker(chunk, key, "whole gene",
                                 domain_motivs, motiv_insertions, 0, 0)

        else:

            sub_domain_cluster, sub_insertions_cluster, domain_part = [], [], []
            subdomain_insert_orient_neg_cluster, subdomain_insert_orient_pos_cluster = 0, 0
            insertions = np.append(insertions, 1)  # avoid overflow

            for i, (domain_insert, sub_domain, sub_insertions, sub_name, neg, pos) \
                in enumerate(zip(insertions, chunk[key].GC_content, chunk[key].motif_seq,
                                 domains, chunk[key].subdomain_insert_orient_neg.values(
                ),
                    chunk[key].subdomain_insert_orient_plus.values())):

                subdomain_insert_orient_neg_cluster += neg
                subdomain_insert_orient_pos_cluster += pos
                sub_domain_cluster.append(sub_domain)
                sub_insertions_cluster.append(sub_insertions)

                # if there is only one domains in the gene
                if sum(sub_domain) == chunk[key].length:
                    chunk = domain_maker(chunk, key, "whole gene", sub_domain, sub_insertions,
                                         subdomain_insert_orient_neg_cluster, subdomain_insert_orient_pos_cluster)
                    break

                # if the domain is the last in the gene.
                if i == len(domains)-1:

                    sub_domain_cluster = np.sum(sub_domain_cluster, axis=0)
                    sub_insertions_cluster = np.sum(
                        sub_insertions_cluster, axis=0)

                    if domains[i]-domain_part[0] == chunk[key].length:
                        domain = "whole gene"

                    else:
                        domain = str(domain_part[0]-start) + \
                            " to " + str(domains[-1]-start)

                    chunk = domain_maker(chunk, key, domain, sub_domain_cluster, sub_insertions_cluster,
                                         subdomain_insert_orient_neg_cluster, subdomain_insert_orient_pos_cluster)
                    subdomain_insert_orient_neg_cluster, subdomain_insert_orient_pos_cluster = 0, 0

                # if the next domain is basically a continuation of the last (also with/out insertions), and it isnt the last domain of the gene
                elif (((domain_insert == 0) & (insertions[i+1] == 0)) or
                      ((domain_insert != 0) & (insertions[i+1] != 0))) and (i < len(domains)-1):
                    domain_part.append(sub_name)

                # if the next domain is different from the previous(domain appending consideration will break), and also not the last domain in the gene
                elif (((domain_insert == 0) & (insertions[i+1] != 0)) or
                      ((domain_insert != 0) & (insertions[i+1] == 0))) and (i < len(domains)-1):

                    domain_part.append(sub_name)

                    sub_domain_cluster = np.sum(sub_domain_cluster, axis=0)
                    sub_insertions_cluster = np.sum(
                        sub_insertions_cluster, axis=0)

                    if domains[i+1]-domain_part[0] == chunk[key].length:
                        domain = "whole gene"
                    else:
                        domain = str(domain_part[0]-start) + \
                            " to " + str(domains[i+1]-start)

                    chunk = domain_maker(chunk, key, domain, sub_domain_cluster, sub_insertions_cluster,
                                         subdomain_insert_orient_neg_cluster, subdomain_insert_orient_pos_cluster)

                    domain_part = [domains[i+1]]
                    sub_domain_cluster, sub_insertions_cluster = [], []
                    subdomain_insert_orient_neg_cluster, subdomain_insert_orient_pos_cluster = 0, 0

                if domain_part[0]-start == domains[-1]-start:
                    break
    return chunk


@njit
def pvaluing_jit(pvaluing_array, pvalue_bol, sig_thresh, i):
    ''' The pvaluing_jit function evaluates the essentiality of genes based on 
    input criteria (p-values, gene length, etc.) and assigns an essentiality 
    score (0-4) to each gene. It also determines whether a gene is too small 
    for assaying, updating the essentiality score accordingly.'''

    gene_half_size = pvaluing_array[i][0] * 0.5
    remove_signal = 0

    if pvalue_bol[i]:  # only appends pvalues that are significant
        if pvaluing_array[i][1] == 1:
            essentiality_score = 4  # essential
        elif pvaluing_array[i][2] >= gene_half_size:
            essentiality_score = 3  # likely essential
        else:
            essentiality_score = 2  # possibly essential
    else:
        essentiality_score = 1  # non-essential

    if (essentiality_score == 1) and (pvaluing_array[i][5] == 0) and (pvaluing_array[i][4] >= sig_thresh):
        remove_signal = 1
        essentiality_score = 0  # too small for assaying

    pvaluing_array[i][3] = essentiality_score

    return remove_signal, pvaluing_array


def pvaluing(names, pvalues_list, pvaluing_array, pvalue, variables, baseline_essentials, baseline_non_essentials):
    ''' The pvaluing function takes gene names, p-values, and other relevant 
    inputs to analyze gene essentiality. It applies FDR correction using the 
    Benjamini-Hochberg procedure to account for multiple testing. 
    The function then updates sets of essential and non-essential genes, 
    as well as rejected baselines for both categories. Finally, it returns 
    the processed p-value array and information about the essential and 
    non-essential genes, along with their respective rejected baselines.'''

    strain_existent_essentials = set()
    strain_existent_nonessentials = set()
    rejected_baseline_essentials = set()
    rejected_baseline_nonessentials = set()

    fdr = statsmodels.stats.multitest.multipletests(
        pvalues_list, pvalue, "fdr_bh")  # multitest correction, fdr_bh

    for i, (name, entry) in enumerate(zip(names, pvaluing_array)):
        remove_signal, pvaluing_array = pvaluing_jit(
            pvaluing_array, fdr[0], fdr[-1], i)

        if remove_signal == 1:
            if (name in baseline_essentials) & (len(baseline_essentials) >= 30):
                baseline_essentials.remove(name)
            elif (name in baseline_non_essentials) & (len(baseline_non_essentials) >= 30):
                baseline_non_essentials.remove(name)

        else:
            if name in baseline_essentials:
                strain_existent_essentials.add(name)
                if fdr[0][i] == False:
                    rejected_baseline_essentials.add(name)

            elif name in baseline_non_essentials:
                strain_existent_nonessentials.add(name)
                if fdr[0][i] == True:
                    rejected_baseline_nonessentials.add(name)

    return pvaluing_array, strain_existent_essentials, fdr[3], rejected_baseline_essentials, \
        rejected_baseline_nonessentials, strain_existent_nonessentials


def basket_storage():
    ''' The basket_storage function forwards the program to the correct annotation
    parser function based on the input annotation file format.'''

    print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.GREEN}INFO{Fore.RESET}] Parsing gene annotations")

    file = variables.annotation_file_paths[0]

    if variables.annotation_type == "gff":
        return gene_info_parser_gff(file)

    elif variables.annotation_type == "gb":
        return gene_info_parser_genbank(file)


def gene_info_parser_genbank(file):
    ''' The gene_info_parser_genbank function takes a genbank file as input and 
    extracts gene information, storing it in a dictionary with Gene class 
    instances as values. It parses the file using the SeqIO module, 
    retrieving attributes such as start, end, orientation, identity, 
    and product for each gene.'''

    basket = {}
    genes = []
    for rec in SeqIO.parse(file, "gb"):
        for feature in rec.features:
            # if feature.type == 'gene':
            start = feature.location.start
            end = feature.location.end
            orientation = feature.location.strand
            
            current = end-start
            if current > variables.biggest_gene:
                variables.biggest_gene = current
            
            try:
                identity = feature.qualifiers['locus_tag'][0]

            except KeyError:
                identity = None

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
                    break  # avoids continuing the iteration and passing to another key, which would make "gene" assume another value
                elif "gene" in key:
                    gene = feature.qualifiers['gene'][0]
                    break  # avoids continuing the iteration and passing to another key, which would make "gene" assume another value
                else:
                    gene = identity

            if identity != None:
                basket[identity] = Gene(gene=gene, start=start, end=end, orientation=orientation,
                                        identity=identity, product=product, contig=rec.id)
                genes.append((start, identity, rec.id))

    basket = inter_gene_annotater(basket, genes)

    for gene in basket:

        domain_size = basket[gene].end - basket[gene].start
        basket[gene].start = int(
            basket[gene].start+(domain_size*variables.subdomain_length[0]))
        basket[gene].end = int(
            basket[gene].start+(domain_size*variables.subdomain_length[1]))

        if basket[gene].orientation == 1:
            basket[gene].orientation = "+"
        else:
            basket[gene].orientation = "-"

    return basket


def inter_gene_annotater(basket, genes):

    genes = list(dict.fromkeys(genes))
    # sort by start position of the gene and contig
    genes.sort(key=lambda x: (x[-1], x[0]))

    count = 0
    for i, gene_entry in enumerate(genes[:-1]):
        contig = gene_entry[-1]
        gene = gene_entry[1]
        gene_next = genes[i+1][1]
        if contig == genes[i+1][-1]:  # same contigs
            gene_down_start_border = variables.ir_size_cutoff + \
                basket[gene].start
            gene_up_start_border = basket[gene_next].end - \
                variables.ir_size_cutoff
            domain_size = gene_up_start_border - gene_down_start_border
            if domain_size >= 1:
                count += 1
                ident = f'IR_{count}_{basket[gene].gene}_{basket[gene].orientation}_UNTIL_{basket[gene_next].gene}_{basket[gene_next].orientation}'
                basket[ident] = Gene(gene=ident,
                                     start=basket[gene].start,
                                     end=basket[gene_next].end,
                                     identity=ident,
                                     contig=contig)

        if contig != genes[i+1][-1]:

            circle_closer = basket[gene].end + variables.ir_size_cutoff
            domain_size = len(variables.genome_seq[contig]) - circle_closer
            if domain_size >= 1:
                count += 1
                ident = f'IR_{count}_{basket[gene].gene}_{basket[gene].orientation}_contig_{contig}_-end'
                basket[ident] = Gene(gene=ident,
                                     start=basket[gene].end,
                                     end=len(variables.genome_seq[contig]),
                                     identity=ident,
                                     contig=contig)

    circle_closer = basket[genes[-1][1]].end + variables.ir_size_cutoff
    domain_size = len(variables.genome_seq[contig]) - circle_closer
    if domain_size > 1:
        count += 1
        ident = f'IR_{count}_{basket[gene].gene}_{basket[gene].orientation}_contig_{contig}_-end'
        basket[ident] = Gene(gene=ident,
                             start=basket[genes[-1][1]].end,
                             end=len(variables.genome_seq[contig]),
                             identity=ident,
                             contig=contig)
        
    domain_size = basket[genes[0][1]].start - variables.ir_size_cutoff
    if domain_size > 1:
        count += 1
        ident = f'IR_{count}_contig_{contig}_-start_{basket[gene].gene}_{basket[genes[0][1]].orientation}'
        basket[ident] = Gene(gene=ident,
                             start=0,
                             end=basket[genes[0][1]].start,
                             identity=ident,
                             contig=contig)
    return basket


def gene_info_parser_gff(file):
    ''' The gene_info_parser_gff function takes a GFF file as input and extracts 
    gene information, storing it in a dictionary with Gene class instances as values. 
    It parses the file line by line, retrieving attributes such as start, end, 
    orientation, identity, and product for each gene with a "CDS" feature.'''

    genes = []
    basket = {}
    with open(file) as current:
        for line in current:
            GB = line.split('\t')  # len(GB)
            if "#" not in GB[0][:3]:  # ignores headers
                if GB[2] == "gene":
                    start = int(GB[3])
                    end = int(GB[4])
                    
                    lenght = end-start
                    if lenght > variables.biggest_gene:
                        variables.biggest_gene = lenght
                    
                    features = GB[8].split(";")  # gene annotation file
                    feature = {}
                    for entry in features:
                        entry = entry.split("=")
                        if len(entry) == 2:
                            feature[entry[0]] = entry[1].replace("\n","")
                    if "gene" in feature:
                        gene = feature["gene"]
                    if "Name" in feature:
                        gene = feature["Name"]
                    else:
                        gene = feature["ID"]
    
                    contig = GB[0]
                    orientation = GB[6]  # orientation of the gene
    
                    basket[feature["ID"]] = Gene(gene=gene, start=start, end=end, orientation=orientation,
                                                 identity=feature["ID"], product=feature["product"], contig=contig)

                    genes.append((start, feature["ID"], contig))

            if "##FASTA" in GB[0]:
                break
    
    if len(basket) == 0:
        print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.YELLOW}WARNING{Fore.RESET}] Watch out, no genomic features were loaded. The gff file is not being parsed correctly.")
        
    basket = inter_gene_annotater(basket, genes)

    for gene in basket:

        domain_size = basket[gene].end - basket[gene].start
        basket[gene].start = int(
            basket[gene].start+(domain_size*variables.subdomain_length[0]))
        basket[gene].end = int(
            basket[gene].start+(domain_size*variables.subdomain_length[1]))

    return basket


def insertions_parser(startup=True):
    ''' This function, insertions_parser, performs the following tasks:

        1. Parses an insertion file to identify unique insertions and their orientation, 
        updating related data structures.

        2. Calculates the total number of insertions and the positive strand transposon 
        insertion ratio.

        3. Compiles transposon motif counts and calculates motif frequencies.

        4. Prints the transposon insertion frequency for each leading strand motif.'''

    # new code
    insertions_df = pd.read_csv(variables.insertion_file_path)
    insertions_df["unique"] = insertions_df["#Contig"] + \
        insertions_df["position"].astype(str)+insertions_df["Orientation"]
    insertions_df.drop_duplicates(subset=['unique'], inplace=True)
    if startup:
        print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.GREEN}INFO{Fore.RESET}] Total Insertions in library: {len(insertions_df)}")
        
    insertions_df.drop(
        ['Read Counts', 'Average mapQ across reads', "unique"], axis=1, inplace=True)

    insertions_df = insertions_df[insertions_df["Relative Position in Gene (0-1)"].astype(
        float) <= variables.subdomain_length[1]]
    insertions_df = insertions_df[insertions_df["Relative Position in Gene (0-1)"].astype(
        float) >= variables.subdomain_length[0]]

    variables.total_insertions = len(insertions_df)

    #borders = insertions_df['Transposon Border Sequence'].tolist()
    orient_pos = len([n for n in insertions_df["Orientation"] if n == "+"])
    variables.positive_strand_tn_ratio = orient_pos / variables.total_insertions

    for (contig, insertion, orientation, border) in zip(insertions_df["#Contig"],
                                                        insertions_df["position"],
                                                        insertions_df["Orientation"],
                                                        insertions_df["Transposon Border Sequence"]):

        if variables.insertions_contig[contig] == {}:
            variables.borders_contig[contig] = [border]
            variables.orientation_contig[contig] = [orientation]
            variables.insertions_contig[contig] = [int(insertion)]

        else:
            variables.borders_contig[contig].append(border)
            variables.orientation_contig[contig].append(orientation)
            variables.insertions_contig[contig].append(int(insertion))
    
    contig_df = insertions_df.groupby("#Contig")
    for (contig,contig_df) in contig_df:
        borders = contig_df['Transposon Border Sequence'].tolist()
        
        variables.transposon_motiv_count[contig] = motiv_compiler(
            borders, variables.regex_compile)
        total = sum(variables.transposon_motiv_count[contig]) #total number of insertions in contig

        variables.transposon_motiv_freq[contig] = np.zeros(
            len(variables.transposon_motiv_count[contig]))
        for i, element in enumerate(variables.transposon_motiv_count[contig]):
            if element != 0:
                variables.transposon_motiv_freq[contig][i] = round(
                    element / total * 100, 1)
    
        # calculating the insertion frequency
        variables.chance_motif_tn[contig] = variables.normalizer(contig)
    
    #if startup:
    #    print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.GREEN}INFO{Fore.RESET}] Transposon insertion frequency (on leading strand):")
    #    for contig in variables.transposon_motiv_freq:
    #        for tn, mtv in zip(variables.transposon_motiv_freq[contig], variables.di_motivs):
    #            print(f" {Fore.GREEN}{mtv}: {Fore.RESET}{tn}%")


def cpu():
    ''' The cpu function determines the available number of CPUs that can be used.
    If more than one is available, one of them is left to avoid clogging the computer,
    while the others are used by the program. '''

    c = multiprocessing.cpu_count()
    if c >= 2:
        c -= 1
    pool = multiprocessing.Pool(processes=c)
    return pool, c


def insertion_annotater(chunk, variables):
    ''' The function insertion_annotater takes a dictionary chunk of gene 
    objects and a set of variables and annotates each gene object with 
    various insertion-related information, such as the number of insertions 
    and their orientation within each subdomain of the gene, the GC content of 
    each subdomain, and the transposon motif content of each subdomain. T
    he function returns the annotated chunk. '''

    for key in chunk:

        chunk[key].subdomain_insert_orient_neg, chunk[key].subdomain_insert_orient_plus = {
        }, {}  # clear past values
        subdomain_insertions, subdomain_insert_seq = {}, {}

        GC_content, domains = [], []
        # first element is start position of gene
        for i, subdomain in enumerate(chunk[key].domains[:-1]):
            GC_content.append(np.array(count_GC([variables.genome_seq[chunk[key].contig][subdomain:chunk[key].domains[i+1]+1]],
                                                variables.regex_compile)))
            domains.append(subdomain)
            # every domain needs an entry, even if empty
            subdomain_insert_seq[subdomain] = [""]
            chunk[key].subdomain_insert_orient_plus[subdomain] = 0
            chunk[key].subdomain_insert_orient_neg[subdomain] = 0

        chunk[key].domain_notes = domains
        chunk[key].GC_content = GC_content
        # pinning all the insertions to their gene domains

        for i, subdomain in enumerate(chunk[key].domains[:-2]):
            domain_start = subdomain-chunk[key].start
            domain_end = chunk[key].domains[i+1]-chunk[key].start
            subdomain_insertions[subdomain] = sum(
                chunk[key].gene_insert_matrix[domain_start:domain_end])
            if subdomain_insertions[subdomain] != 0:
                for j in range(subdomain-1, chunk[key].domains[i+1]-1):
                    border = variables.borders_contig[chunk[key].contig][j]
                    insert = variables.insertions_contig[chunk[key].contig][j]
                    if insert == 1:
                        subdomain_insert_seq[subdomain] = subdomain_insert_seq[subdomain] + [
                            border]
                        # positive insertion
                        if variables.orientation_contig_plus[chunk[key].contig][j] == 1:
                            chunk[key].subdomain_insert_orient_plus[subdomain] += 1
                        # negative insertion
                        if variables.orientation_contig_neg[chunk[key].contig][j] == 1:
                            chunk[key].subdomain_insert_orient_neg[subdomain] += 1

        # Transposon motif content in each domain (alwyas 16, corresponding to the dinucleotide combo)
        if subdomain_insert_seq != {}:
            motif_seq = []
            for key1, value in sorted(subdomain_insert_seq.items()):
                motif_seq.append(motiv_compiler(
                    value, variables.regex_compile))
            chunk[key].motif_seq = motif_seq

        else:
            chunk[key].motif_seq = np.zeros(16)

        chunk[key].domain_insertions_total = np.zeros(
            len(subdomain_insertions))
        for i, (key1, value) in enumerate(sorted(subdomain_insertions.items())):
            chunk[key].domain_insertions_total[i] = value

    chunk = essentials(chunk, variables)

    return chunk


def multi_annotater(basket):
    ''' The function insertion_annotater takes a dictionary chunk of gene 
    objects and a set of variables and annotates each gene object with 
    various insertion-related information, such as the number of insertions 
    and their orientation within each subdomain of the gene, the GC content of 
    each subdomain, and the transposon motif content of each subdomain. T
    he function returns the annotated chunk. '''

    pool, cpus = cpu()
    # divides the dictionary keys into smaller blocks that can be efficiently be multiprocessed
    divider = len(basket)//cpus
    return_list = [dict() for i in range(cpus)]
    i, list_iter = 0, 0
    for k in basket:
        if i < cpus:
            return_list[i][k] = basket[k]

        if i == cpus:  # odd number split will be distributed equally
            list_iter += 1
            return_list[list_iter][k] = basket[k]

        elif len(return_list[i]) >= divider:
            i += 1

    result_objs = []
    for chunk in return_list:
        result = pool.apply_async(
            insertion_annotater, args=((chunk, variables)))
        result_objs.append(result)

    pool.close()
    pool.join()
    result = [result.get() for result in result_objs]
    # demultiplexing the results
    for subresult in result:
        for key in basket:
            if key in subresult:
                basket[key] = subresult[key]

    return basket


def pvalue_iteration(names, pvalues_list, pvaluing_array, pvalue, pvalue_listing, euclidean_points, variables):
    ''' This function performs a p-value iteration analysis to determine the 
    true positive rate (TPR) and specificity of a given test at different 
    p-value thresholds. It calls the function 'pvaluing' to calculate the 
    number of true positives, false negatives, true negatives and false 
    positives at a given p-value threshold. It then calculates the TPR and 
    specificity from these values and appends them, along with the p-value 
    threshold, to two lists. Finally, it returns these two lists. '''

    pvaluing_array_copy = pvaluing_array.copy()
    baseline_essentials = set(variables.true_positives.keys())
    baseline_non_essentials = set(variables.true_negatives.keys())
    
    pvaluing_array_discard, strain_existent_essentials, fdr, rejected_baseline_essentials, \
        rejected_baseline_nonessentials, strain_existent_nonessentials = pvaluing(
            names, pvalues_list, pvaluing_array_copy, pvalue, variables, baseline_essentials, baseline_non_essentials)

    if len(pvaluing_array_discard) == 1:  # if there is only the header
        rejected_baseline_essentials = strain_existent_essentials

    TP = len(strain_existent_essentials) - len(rejected_baseline_essentials)
    FN = len(rejected_baseline_essentials)

    if (TP + FN) != 0:
        TPR_sensitivity = TP / float(TP + FN)
    else:
        TPR_sensitivity = 0

    TN = len(strain_existent_nonessentials) - \
        len(rejected_baseline_nonessentials)
    FP = len(rejected_baseline_nonessentials)

    if (TN + FP) != 0:
        specificity = 1 - (TN / float(TN + FP))
    else:
        specificity = 0

    pvalue_listing.append(pvalue)
    euclidean_points.append([specificity] + [TPR_sensitivity])

    return pvalue_listing, euclidean_points


def class_to_numba(basket):
    ''' The class_to_numba function takes a dictionary basket containing gene 
    information and extracts significant gene domains to convert them into 
    arrays that can be used with the numba library. It first loops through the 
    dictionary to extract gene names, p-values, and domain-related information. 
    Then, it creates a pvaluing_array containing this information in a numerical 
    format. Finally, it calls the reference_essentials_loader function to load 
    essential and non-essential gene lists, which are returned along with the 
    pvaluing_array, names, and pvalues_list.  '''

    pvaluing_array, names, pvalues_list = [], [], []
    for gene in basket:
        for domain in basket[gene].significant:
            part = 0
            if basket[gene].significant[domain].domain_part == "whole gene":
                part = 1
            names.append(basket[gene].gene)
            pvalues_list.append(basket[gene].significant[domain].pvalue)
            pvaluing_array.append(np.array((basket[gene].length,
                                           part,
                                           basket[gene].significant[domain].dom_len,
                                           1,  # non essential
                                           basket[gene].significant[domain].pvalue,
                                           basket[gene].significant[domain].dom_insert),
                                           dtype=np.float64))

    pvaluing_array = np.array(pvaluing_array)
    return pvaluing_array, names, pvalues_list


def multi_pvalue_iter(basket):
    ''' This function performs multiple iterations of p-value thresholding to 
    determine true positive and false positive rates for a set of genes, 
    using parallel processing for efficiency. It takes a dictionary of genes 
    and their properties, converts them to Numpy arrays, and returns the 
    resulting lists of p-value thresholds and corresponding true/false positive rates.  '''

    # change here to increase resolution of iteration
    pvalue = [variables.pvalue * 0.5 ** i for i in range(200)]
    pool, cpus = cpu()
    result_objs, pvalue_listing, euclidean_points = [], [], []

    pvaluing_array, names, pvalues_list = class_to_numba(basket)

    for p in pvalue:
        result = pool.apply_async(pvalue_iteration, args=((names, pvalues_list, 
                                                           pvaluing_array, p, pvalue_listing,
                                                           euclidean_points, variables)))
        result_objs.append(result)

    pool.close()
    pool.join()

    result = [result.get() for result in result_objs]

    result = sorted(result, key=lambda e: e[0][0], reverse=True)
    pvalue_listing, euclidean_points = zip(*result)

    return pvalue_listing, [n for n in zip(*euclidean_points)]


def final_compiler(optimal_basket, pvalue, euclidean_points):
    ''' The final_compiler function takes a dictionary of genes and their 
    properties and performs essentiality calling and statistical analysis 
    using p-value thresholding. The function converts the gene dictionary into 
    Numpy arrays and applies the pvaluing function to determine essentiality. 
    The function then compiles a table of gene properties and essentiality 
    information, and generates an ROC (Receiver Operating Characteristic) 
    curve plot for the data. Finally, the function saves the table and plot 
    to files in a specified directory. '''

    pvaluing_array, names, pvalues_list = class_to_numba(optimal_basket)

    fdr = statsmodels.stats.multitest.multipletests(
        pvalues_list, pvalue, "fdr_bh")  # multitest correction, fdr_bh

    for i, (name, entry) in enumerate(zip(names, pvaluing_array)):
        remove_signal, pvaluing_array = pvaluing_jit(
            pvaluing_array, fdr[0], fdr[-1], i)
    
    fdr = fdr[3]
    legenda = f"Threshold p-value: {pvalue}"
    essentiality = []

    for i, entry in enumerate(pvaluing_array):
        if entry[3] == 0:
            essentiality.append("too small for assaying")
        if entry[3] == 1:
            essentiality.append("Non-Essential")
        if entry[3] == 2:
            essentiality.append("Probably Essential")
        if entry[3] == 3:
            essentiality.append("Likelly Essential")
        if entry[3] == 4:
            essentiality.append("Essential")

    genes_list = [["Total unique Tn insertions"] + ["Essentiality p-value"] +
                  ["Contig"] + ["Gene ID"] + ["Description"] + ["Domain length"] +
                  ["Insertion Orientation ratio (+/total)"] + ["Insertion Orientation p-value"] +
                  ["gene region"] + ["Gene Orientation"] + ["Gene name"] + ["Essentiality"]]

    i = 0
    for gene in optimal_basket:
        for domain in optimal_basket[gene].significant:
            genes_list.append([optimal_basket[gene].significant[domain].dom_insert] +
                              [str(optimal_basket[gene].significant[domain].pvalue).replace("'","")] +
                              [optimal_basket[gene].contig] +
                              [optimal_basket[gene].identity] +
                              [optimal_basket[gene].product] +
                              [int(optimal_basket[gene].significant[domain].dom_len/2)] +
                              [optimal_basket[gene].significant[domain].ratio_orient] +
                              [optimal_basket[gene].significant[domain].orient_pvalue] +
                              [optimal_basket[gene].significant[domain].domain_part] +
                              [optimal_basket[gene].orientation] +
                              [optimal_basket[gene].gene] +
                              [essentiality[i]])
            i += 1

    # there is some issue with the couting
    significant_genes_list = list(
        filter(lambda x: x[-1] == "Essential", genes_list))  # gets just the essentials
    significant_genes_list_full = significant_genes_list + \
        list(filter(lambda x: x[-1] == "Likelly Essential",
             genes_list))  # gets just the essentials
    significant_genes_list_full = significant_genes_list_full + \
        list(filter(lambda x: x[-1] == "Possibly Essential",
             genes_list))  # gets just the essentials
    # gets just the non evaluated genes
    non_assayed = list(
        filter(lambda x: x[-1] == "too small for assaying", genes_list))
    # gets just the non essentials
    non_essentials = list(
        filter(lambda x: x[-1] == "Non-Essential", genes_list))

    essentials = set()
    a = [essentials.add(gene[3]) for gene in significant_genes_list_full]

    non_assayed_list = set()
    a = [non_assayed_list.add(gene[3])
         for gene in non_assayed if gene[-2] not in non_assayed_list]

    non_essentials_list = set()
    a = [non_essentials_list.add(
        gene[3]) for gene in non_essentials if gene[-2] not in non_essentials_list]

    intersect = non_assayed_list.intersection(essentials)
    intersect |= non_assayed_list.intersection(non_essentials_list)
    full_na_genes = len(non_assayed_list) - len(intersect)
    intersect = non_essentials_list.intersection(essentials)
    intersect |= non_essentials_list.intersection(non_assayed_list)
    full_non_e_genes = len(non_essentials_list) - len(intersect)

    dic = {}
    for contig in variables.transposon_motiv_freq:
        dic[contig] = {}
        for x, i in zip(variables.transposon_motiv_freq[contig], variables.di_motivs):
            dic[contig][i] = x

    genes_list.insert(
        0, ["#Transposon insertion percent bias (+strand): %s" % dic])
    genes_list.insert(0, ["#p-value cutoff: %s" % pvalue])
    genes_list.insert(0, ["#fdr corrected p-value cutoff: %s" % fdr])
    genes_list.insert(
        0, ["#Number of genes with at least one domain that is non-essential: %s" % len(non_essentials_list)])
    genes_list.insert(
        0, ["#Number of whole genes that are non-essential: %s" % full_non_e_genes])
    genes_list.insert(
        0, ["#Number of genes with at least one domain too small for assaying: %s" % len(non_assayed_list)])
    genes_list.insert(
        0, ["#Number of whole genes too small for assaying: %s" % full_na_genes])
    genes_list.insert(
        0, ["#Number of genes with at least one domain that is essential: %s" % len(essentials)])
    genes_list.insert(
        0, ["#Number of whole genes that are essential: %s" % len(significant_genes_list)])

    output_writer(variables.directory, variables.output_name, genes_list)
    x, y = zip(*euclidean_points[0])

    fig, ax1 = plt.subplots()
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    ax1.set_xlabel("True Negative Rate")  # 1 - specificity
    ax1.set_ylabel("True Positive Rate")  # (sensitivity)
    ax1.plot(x, y, color="midnightblue")
    ax1.plot(ax1.get_xlim(), ax1.get_ylim(), ls="--", c=".3")
    ax1.tick_params(axis='y')
    ax1.set_title("Reciver Operator Curve (ROC)\nused to auto-determine the essentiality calling threshold for %s" %
                  variables.strain, size=9, pad=13)
    ax1.legend([legenda], loc="lower right")
    fig.tight_layout()
    plt.savefig(
        f"{variables.directory}/ROC_curves{variables.strain}.png", dpi=300)
    plt.close()
    gene_essentiality_compressor(optimal_basket,genes_list)


def genome_loader(startup=True):
    ''' This function loads the genome sequence and stores it in a dictionary 
    named `genome_seq` under their respective contigs. It also initializes 
    empty dictionaries named `borders_contig`, `orientation_contig`, and 
    `insertions_contig` for each contig. The function works differently 
    depending on the format of the annotation file: if it is in GFF format, 
    it reads the sequence from a FASTA file referenced in `annotation_file_paths[1]`,
    while if it is in GenBank format, it reads the sequence from the GenBank 
    file referenced in `annotation_file_paths[0]`. The function also prints 
    out the name of each loaded contig. '''

    if variables.annotation_type == "gff":
        with open(variables.annotation_file_paths[1]) as current:
            for line in current:
                if ">" in line:
                    contig = line.split()[0][1:]
                    variables.genome_seq[contig] = ""
                    variables.borders_contig[contig] = {}
                    variables.orientation_contig[contig] = {}
                    variables.insertions_contig[contig] = {}
                else:
                    variables.genome_seq[contig] += line[:-1]

        for contig in variables.genome_seq:
            if startup:
                print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.GREEN}INFO{Fore.RESET}] Loaded contig {contig}")
            variables.genome_length += len(variables.genome_seq[contig])

    elif variables.annotation_type == "gb":
        for rec in SeqIO.parse(variables.annotation_file_paths[0], "gb"):
            variables.genome_length = len(rec.seq)
            variables.annotation_contig = rec.id
            variables.borders_contig[variables.annotation_contig] = {}
            variables.orientation_contig[variables.annotation_contig] = {}
            variables.insertions_contig[variables.annotation_contig] = {}
            if startup:
                print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.GREEN}INFO{Fore.RESET}] Loaded contig {variables.annotation_contig}")
            variables.genome_seq[variables.annotation_contig] = str(rec.seq)


def domain_iterator(basket):
    ''' The function domain_iterator performs an iterative process to find 
    the optimal size of genomic regions to divide a genome for downstream analysis. 
    It starts by defining a maximum size of the genomic region, then iteratively
    resizes the regions and performs a ROC (receiver operating characteristic) 
    analysis to determine the optimal p-value threshold for the current region 
    size. The function uses the multi_annotater and domain_resizer functions to 
    prepare the genomic regions for the ROC analysis. The output is the optimized 
    set of genomic regions, the optimal p-value threshold, and the ROC data used 
    to determine the optimal threshold. '''

    def ROC(basket):
        def orientation_pvaluing(basket):
            pvalues_orient = []
            for key in basket:
                for domain in basket[key].significant:
                    pvalues_orient.append(
                        basket[key].significant[domain].orient_pvalue)
            pvalues_orient = [
                score for score in pvalues_orient if score != "N/A"]

            if len(pvalues_orient) > 0:
                fdr = statsmodels.stats.multitest.multipletests(
                    pvalues_orient, 0.01, "fdr_bh")  # multitest correction, fdr_bh
            p = 0
            for key in basket:
                for domain in basket[key].significant:
                    if basket[key].significant[domain].orient_pvalue != "N/A":
                        basket[key].significant[domain].orient_pvalue = fdr[1][p]
                        p += 1
            return basket

        basket = orientation_pvaluing(basket)
        pvalue_listing, euclidean_points = multi_pvalue_iter(basket)
        # point where TPR_sensitivity = 1 and 1-specificity = 0
        convergence = [[0] + [1]] * len(np.array(euclidean_points[0]))

        #lower is better
        distances = scipy.spatial.distance.cdist(
            np.array(euclidean_points[0]), np.array(convergence), 'euclidean')

        optimal_distance = []
        # returns first entry of the array of ecludian
        list(filter(lambda x: optimal_distance.append(x[0]), distances))
        # gives the index of the inflexion point. that minimizes the eucledian distances
        inflexion_points_index = (
            optimal_distance.index(min(optimal_distance)))

        pvalue = pvalue_listing[inflexion_points_index][0]
        best_distance = optimal_distance[inflexion_points_index]

        return pvalue, best_distance, np.array(euclidean_points)

    def iterating(i, iterator_store, euclidean_distances, basket, current_gap):
        basket = domain_resizer(variables.domain_iteration[i], basket)
        basket = multi_annotater(basket)
        best_pvalue, best_distance, euclidean_points = ROC(basket)
        iterator_store.append(
            [current_gap] + [best_pvalue] + [best_distance] + [i] + [euclidean_points])
        return iterator_store

    euclidean_distances = []
    iterator_store = []
    i, current_gap = 0, 0

    try:

        iteration = 1
        for f in variables.domain_iteration:
            if int(variables.genome_length / variables.total_insertions * f) <= variables.biggest_gene:
                iteration += 1
    
        while (current_gap <= variables.biggest_gene) and (i+1 < len(variables.domain_iteration)):
            current_gap = int(variables.genome_length /
                              variables.total_insertions * variables.domain_iteration[i])
            print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.GREEN}INFO{Fore.RESET}] Current domain division iteration: {i+1} out {iteration}")
            iterator_store = iterating(
                i, iterator_store, euclidean_distances, basket, current_gap)
            i += 1
            
        # sort by best eucledian distance
        sorted_optimal = sorted(iterator_store, key=lambda e: e[2])
        variables.best_domain_size = sorted_optimal[0][0]
        
    except Exception:
        print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.YELLOW}WARNING{Fore.RESET}] Low saturating library detected")
        variables.best_domain_size = sorted_optimal[0]
    
    print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{Fore.GREEN}INFO{Fore.RESET}] Optimal domain division size of {variables.best_domain_size}bp")
    fig, ax1 = plt.subplots()

    plt.xlim(0, 1)
    plt.ylim(0, 1)

    ax1.set_xlabel("True Negative Rate")  # 1 - specificity
    ax1.set_ylabel("True Positive Rate")  # (sensitivity)

    legenda = [n[0] for n in iterator_store]
    for i in iterator_store:
        z, x = zip(*i[4][0])
        ax1.plot(z, x)

    ax1.plot(ax1.get_xlim(), ax1.get_ylim(), ls="--", c=".3")
    ax1.tick_params(axis='y')
    ax1.set_title("Reciver Operator Curve (ROC)\nused to auto-determine the essentiality calling threshold for %s" %
                  variables.strain, size=9, pad=13)
    ax1.legend(legenda, loc="lower right")
    fig.tight_layout()
    plt.savefig(os.path.join(variables.directory,
                "ROC_curves_iterator_%s.png" % variables.strain), dpi=300)
    plt.close()
    best_index = int(sorted_optimal[0][3])
    output_writer(variables.directory, "Best_ROC_points{}".format(
        variables.output_name), iterator_store[best_index][4][0])

    return best_index, iterator_store[best_index][1], iterator_store[best_index][4]

def gene_essentiality_compressor(basket,genes_list):
    
    essentials = pd.read_csv(f"{variables.directory}/{variables.output_name}.csv",
                             skiprows=9,low_memory=False)

    classifier = {'Essential': 1,
                  'Likelly Essential': 2,
                  'Probably Essential': 3,
                  'too small for assaying': 4,
                  'Non-Essential': 5}

    binned_essentiality = {}
    binned_essentiality_domain = {}

    for gene_name, essential, domain in zip(essentials['Gene ID'], essentials['Essentiality'], essentials['gene region']):
        if essential in classifier:
            essential = classifier[essential]

        size = 0
        if domain != 'whole gene':
            initial, end = domain.split(" to ")
            size = int(end) - int(initial)

        if gene_name in binned_essentiality:
            binned_essentiality[gene_name].add(essential)

        else:
            binned_essentiality[gene_name] = set([essential])
            binned_essentiality_domain[gene_name] = {essential: 0}
        # to sort if a gene has more domains considered non essential or not
        if essential in binned_essentiality_domain[gene_name]:
            binned_essentiality_domain[gene_name][essential] += size
        else:
            binned_essentiality_domain[gene_name] = {
                **binned_essentiality_domain[gene_name], **{essential: size}}

    for key in binned_essentiality_domain:
        if (4 in binned_essentiality_domain[key]) and \
            (5 in binned_essentiality_domain[key]) and \
                (len(binned_essentiality_domain[key]) == 2):
            # if there are more domains that are too small than the sum of all domains
            if binned_essentiality_domain[key][4] > (binned_essentiality_domain[key][5] + binned_essentiality_domain[key][4]) * variables.domain_uncertain_threshold:
                binned_essentiality[key] = {4}
            else:
                binned_essentiality[key] = {5}

    essentials['Full gene essentiality'] = None
    for gene_name, i in zip(essentials['Gene ID'], essentials.index):
        if gene_name in binned_essentiality:
            for j in range(1, 6):  # from the classifier dict
                if j in binned_essentiality[gene_name]:
                    essentials.loc[i, 'Full gene essentiality'] = [
                        n for n in classifier if classifier[n] == j][0]
                    break
    
    final_out = [["Gene ID"] + ["Gene Name"] + ["Gene Orientation"] +
                  ["Contig"] + ["Gene Lenght"] + ["Gene Product"] +
                  ["Total Number of insertions"] + ["Essentiality"]]

    for gene in basket:
        final_out.append([basket[gene].identity] +
                          [basket[gene].gene] +
                          [basket[gene].orientation] +
                          [basket[gene].contig] +
                          [basket[gene].length] +
                          [basket[gene].product] +
                          [sum(basket[gene].domain_insertions_total)] +
                          [essentials[essentials["Gene ID"]==gene]["Full gene essentiality"].iloc[0]])
    
    output_writer(variables.directory, f"{variables.output_name}_final", final_out)

def main(argv):
    ''' This function is the main function that calls all the other 
    functions in the program. It takes in a list of command-line arguments as 
    input (argv), sets up the necessary inputs and paths, loads the genome 
    sequence data, parses the transposon insertion data, constructs a basket 
    of genes, generates a matrix of gene insertions, iterates over different 
    domain sizes and identifies the optimal one, compiles the final results, 
    and prints the start and end times of the program execution. '''
    
    inputs(argv)
    path_finder()
    genome_loader(startup=True)
    insertions_parser(startup=True)
    evaluator_basket = blast_maker()
    if (len(variables.true_positives) != 0) or (len(variables.true_negatives) != 0):
        evaluator_basket = gene_insertion_matrix(evaluator_basket)
        best_index, pvalue, euclidean_points = domain_iterator(evaluator_basket)
    else:
        best_index, pvalue, euclidean_points = 0,variables.pvalue,[[[0,0],[0,0]]]

    genome_loader(startup=False) # create a cleaner version where this doesnt need to repeat
    insertions_parser(startup=False)
    basket = basket_storage()
    basket = gene_insertion_matrix(basket)
    best_basket = domain_resizer(variables.domain_iteration[best_index], basket)
    best_basket = multi_annotater(best_basket)
    final_compiler(best_basket, pvalue, euclidean_points)


if __name__ == "__main__":
    if len(sys.argv) > 2:
        argv = sys.argv[1:]
    main(argv)
    multiprocessing.set_start_method("spawn")
