import statsmodels.stats.multitest
from Bio import SeqIO
import os
import numpy as np
import scipy
import re
from tnseeker.extras.possion_binom import PoiBin
from tnseeker.extras.helper_functions import colourful_errors,csv_writer,subprocess_cmd,file_finder,variables_parser,gb_parser,gff_parser,inter_gene_annotater
from scipy.stats import binomtest
import multiprocessing
import matplotlib.pyplot as plt
import sys
from numba import njit
import importlib.resources as resources
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Dict, List
import matplotlib
matplotlib.use('Agg')

"""
Disclaimer: 
    I am aware this script is messy and complex. 
    It was 6y in the making and changed a lot over time to accommodate endless 
    features, exceptions, and use cases. 
    I learnt Python by building this script. AND IT SHOWS!
    It's more of a pathwork than anything else at the moment.
    Please don't judge me too harshly.
"""

def inputs(argv):
    ''' The inputs function initializes a global Variables instance with various attributes 
    from a given list of command-line arguments. It sets up parameters for directory, 
    strain, annotation, subdomain length, p-value, and output, as well as true positives 
    and true negatives file paths. '''

    global variables
    variables = Variables()
    variables_input = variables_parser(argv)
    variables.directory = variables_input["directory"]
    variables.strain = variables_input["strain"]
    variables.annotation_type = variables_input["annotation_type"]
    variables.annotation_folder = variables_input["annotation_folder"]
    variables.subdomain_length = [float(variables_input["subdomain_length_up"]), float(variables_input["subdomain_length_down"])]
    variables.pvalue = float(variables_input["pvalue"])
    variables.ir_size_cutoff = int(variables_input["intergenic_size_cutoff"])
    variables.output_name = variables.strain + "_essential_features"
    variables.domain_uncertain_threshold = float(variables_input["domain_uncertain_threshold"])
    variables.biggest_gene = 0

    variables.cpus = int(variables_input["cpus"])

    variables.true_positives = resources.files('tnseeker').joinpath('data/true_positives.fasta')
    variables.true_negatives = resources.files('tnseeker').joinpath('data/true_negatives.fasta')


def path_finder():
    ''' The path_finder function searches for and sets file paths for both insertions file and the 
    correct annotation (gb/gff) based on the global variables instance. It uses the nested sub_path_finder 
    function to search for specific file types and names within specified folders. '''

    variables.insertion_file_path = file_finder(variables.directory, ['*.csv'], "all_insertions")

    if variables.annotation_type == "gb":
        extention = ['*.gb']
        if file_finder(variables.annotation_folder, extention, variables.strain) == None:
            extention = ['*.gbk']
        variables.annotation_file_paths = [file_finder(
            variables.annotation_folder, extention, variables.strain)]

    elif variables.annotation_type == "gff":
        variables.annotation_file_paths = [file_finder(
            variables.annotation_folder, ['*.gff'], variables.strain)]
        variables.annotation_file_paths.append(file_finder(
            variables.annotation_folder, ['*.fasta'], variables.strain))

        if len(variables.annotation_file_paths) < 2:
            colourful_errors("FATAL",
                 "Check there is a .fasta and .gff file in the indicated folder.")
            raise Exception

def output_writer(output_folder, name_folder, output_file):
    output_file_path = os.path.join(output_folder, name_folder + ".csv")
    csv_writer(output_file_path,output_file)

@dataclass
class Variables:
    """ 
    A comprehensive container for various attributes and methods related to genomic data analysis.
    """

    directory: Optional[str] = None
    strain: Optional[str] = None
    annotation_type: Optional[str] = None
    annotation_folder: Optional[str] = None
    pan_annotation: Optional[str] = None
    output_name: Optional[str] = None
    true_positives: Optional[int] = None
    true_negatives: Optional[int] = None
    insertion_file_path: Optional[str] = None
    annotation_file_paths: Optional[List[str]] = None
    borders_contig: Dict = field(default_factory=dict)
    orientation_contig: Dict = field(default_factory=dict)
    insertions_contig: Dict = field(default_factory=dict)
    genome_seq: Dict = field(default_factory=dict)
    genome_length: int = 0
    annotation_contig: Optional[Dict] = None
    total_insertions: Optional[int] = None
    positive_strand_tn_ratio: Optional[float] = None
    transposon_motiv_count: Dict = field(default_factory=dict)
    transposon_motiv_freq: Dict = field(default_factory=dict)
    chance_motif_tn: Dict = field(default_factory=dict)
    orientation_contig_plus: Dict = field(default_factory=dict)
    orientation_contig_neg: Dict = field(default_factory=dict)
    subdomain_length: List[float] = field(default_factory=lambda: [0.1, 0.9])
    pvalue: float = 0.1
    domain_iteration:  List[int] = field(init=False, repr=False)
    regex_compile: List[re.Pattern] = field(init=False, repr=False)
    di_motivs: List[str] = field(init=False, repr=False)

    def __post_init__(self):
        self.regex_compile, self.di_motivs = self.motif_compiler()
        self.domain_iteration = self.fibonacci(1, 1, 50000, [])

    def fibonacci(self,a,b,n,result):
        c = a+b
        result.append(c)
        if c < n:
            self.fibonacci(b,c,n,result)
        return result

    def normalizer(self,contig):
        motif_genome = np.zeros(16)
        # Counts the entire motifv content of the genome; determining the probability of having an insertion in each motif
        atcg = count_GC([self.genome_seq[contig]], self.regex_compile)
        for i, element in enumerate(motif_genome):
            motif_genome[i] = element + atcg[i]

        # probability of having a motif with a transposon in it
        self.chance_motif_tn[contig] = np.divide(
            self.transposon_motiv_count[contig], motif_genome)
        self.chance_motif_tn[contig] = np.clip(self.chance_motif_tn[contig], 0, 0.999999)
        return np.array([[n] for n in self.chance_motif_tn[contig]])

    def motif_compiler(self):
        dna = ["A", "T", "C", "G"]
        di_motivs = [f"{a}{b}" for a in dna for b in dna]
        prog = [re.compile(f"(?=({m}))", re.IGNORECASE) for m in di_motivs]
        return prog, di_motivs

@dataclass
class Significant:
    """ 
    A container for storing information about a gene or domain.
    
    Attributes:
        dom_insert (Optional[float]): Location matrix of domain insertions.
        pvalue (Optional[float]): computed essentiality p-value.
        dom_len (Optional[int]): Domain length.
        ratio_orient (Optional[float]): matrix of orientation ratio.
        orient_pvalue (Optional[float]): Orientation p-value.
        domain_part (Optional[str]): Gene region.
        essentiality (str): Essentiality classification (default: "Non-Essential").
    
    Used for the final processing of essentiality.
    """

    dom_insert: Optional[float] = None
    pvalue: Optional[float] = None
    dom_len: Optional[int] = None
    ratio_orient: Optional[float] = None
    orient_pvalue: Optional[float] = None
    domain_part: Optional[str] = None
    essentiality: str = "Non-Essential"

@dataclass
class Gene:
    """ 
    A container for storing information about a gene, including its attributes such as:
    - Start and end positions
    - Orientation, domains, and identity
    - Gene product and associated contig
    - Various data related to insertions, GC content, and domain significance
    
    This class is used as the first step for storing processed gene information from annotation files.
    """

    gene: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    start_trim: Optional[int] = None
    end_trim: Optional[int] = None
    orientation: Optional[str] = None
    domains: Optional[str] = None
    identity: Optional[str] = None
    product: Optional[str] = None
    contig: Optional[str] = None
    gene_insert_matrix: Optional[str] = None  # Assuming matrix is stored as a string or another type
    domain_insertions_total: Optional[int] = None
    GC_content: Optional[float] = None
    domain_notes: Optional[str] = None
    motif_seq: Optional[str] = None
    subdomain_insert_orient_plus: Dict = field(default_factory=dict)
    subdomain_insert_orient_neg: Dict = field(default_factory=dict)
    significant: Dict = field(default_factory=dict)
    length: Optional[int] = field(init=False, repr=False)

    def __post_init__(self):
        self.length = (self.end - self.start) if self.start is not None and self.end is not None else None

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
                                           start_trim=start,
                                           end_trim=end,
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
    
    colourful_errors("INFO",
        f"Found {len(variables.true_positives)} benchmark essential genes.")

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
    
    if domain_size == 0:
        domain_size = 1

    for key in basket:
        local_stop = [basket[key].start_trim]
        start_iterator = basket[key].start_trim
        end = basket[key].end_trim

        divider = int((end - start_iterator) / domain_size)

        if divider >= 1:
            for i in range(divider):
                local_stop.append(start_iterator+domain_size)
                start_iterator += domain_size

        if (end - (local_stop[-1] + domain_size)) > domain_size:
            del local_stop[-1]

        if local_stop[-1] != end:
            local_stop.append(end)
        # the end of the gene is always required in duplicate due to poorly optimized downstream functions
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

    colourful_errors("INFO",
        "Compiling insertion matrix.")

    for contig in variables.genome_seq:
        contig_size = len(variables.genome_seq[contig])
        genome_insert_matrix = np.zeros(contig_size, dtype=np.int16)
        genome_borders_matrix = np.zeros(contig_size, dtype='<U2')
        genome_orient_plus_matrix = np.zeros(contig_size, dtype=np.int16)
        genome_orient_neg_matrix = np.zeros(contig_size, dtype=np.int16)
        variables.orientation_contig_plus[contig] = {}
        variables.orientation_contig_neg[contig] = {}

        if len(variables.insertions_contig[contig]) != 0:
            inserts = np.array(variables.insertions_contig[contig], dtype=int)
            borders = np.array(variables.borders_contig[contig])
            orient = np.array(variables.orientation_contig[contig])
        else:
            inserts = np.zeros(contig_size, dtype=np.int16)
            borders = np.zeros(contig_size, dtype='<U2')
            orient = np.zeros(contig_size, dtype=np.int16)

        variables.insertions_contig[contig],\
            variables.borders_contig[contig],\
            variables.orientation_contig_plus[contig],\
            variables.orientation_contig_neg[contig] = \
            contig_matrix_generator(orient, genome_orient_plus_matrix, genome_orient_neg_matrix,
                                    borders, genome_borders_matrix, inserts, genome_insert_matrix, contig)

    for key in basket:
        start = basket[key].start_trim
        end = basket[key].end_trim
        for contig in variables.insertions_contig:
            # calculating local insertion transposon density in bp windows (size variable)
            if basket[key].contig == contig:
                basket[key].gene_insert_matrix = variables.insertions_contig[contig][start-1:end-1]

    return basket


def count_GC(seq, prog):
    """
    Count the occurrences of each regex pattern across a list of sequences.

    Parameters
    ----------
    seq : list of str
        A list of DNA/RNA or other text-based sequences to be searched.
    prog : list of re.Pattern
        A list of compiled regex patterns.

    Returns
    -------
    list of int
        A list of counts. For each pattern in `prog`, and for each string in `seq`,
        the function appends the total number of matches found.
    """
    result = []
    for motiv in prog:
        for insert in seq:
            result.append(len([m.start() for m in re.finditer(motiv, insert)]))
    return result


def motiv_compiler(seq, prog):
    """
    Count how many entries in a list of sequences match each regex pattern, 
    stopping at the first match per sequence.

    Parameters
    ----------
    seq : list of str
        A list of sequences to be checked.
    prog : list of re.Pattern
        A list of compiled regex patterns.

    Returns
    -------
    numpy.ndarray
        A 1D array of length 16 (by default) where each index corresponds to a pattern 
        in `prog`. The value at each index represents how many sequences matched 
        that pattern first. Only the first matching pattern is counted per sequence.
    """
    motiv_inbox = np.zeros(16)
    for insertion in seq:
        if insertion != "":
            for i, motiv in enumerate(prog):
                if [m.start() for m in re.finditer(motiv, insertion)] != []:
                    motiv_inbox[i] = motiv_inbox[i] + 1
                    break
    return motiv_inbox


def poisson_binomial(events, motiv_inbox, tn_chance):
    """
    Compute the Poisson binomial cumulative distribution function (CDF) for 
    a set of events, motif occurrences, and transposon insertion probabilities.

    Parameters
    ----------
    events : numpy.ndarray
        An array indicating how many times each event occurs.
    motiv_inbox : numpy.ndarray
        An array counting motif occurrences (or successful matches).
    tn_chance : numpy.ndarray
        An array of probabilities (chances) for transposon insertion events.

    Returns
    -------
    float
        The absolute value of the Poisson binomial CDF evaluated at the 
        total number of successes, where 'successes' is the sum of motif 
        transposon insertion occurrences (capped by the sum of events).

    Notes
    -----
    - The function defines an inner helper, `chance_matrix`, which constructs a 
    probability array (`p_effective`) for each event occurrence.
    - `PoiBin` is used to calculate the Poisson binomial distribution.
    - The return value is wrapped in `abs()` to handle potential negative 
    rounding issues in extreme cases.
    """
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
            colourful_errors("WARNING",
                f"Due to {e}, {key} could not be evaluated.")
            chunk[key].significant[domain].pvalue = 1
        return chunk

    # annexes all the genes which fit the criteria of essentiality
    for key in chunk:

        chunk[key].significant = {}  # clears any overlapping dictionaries
        # total number of insertions in gene
        insertions = chunk[key].domain_insertions_total
        domains = chunk[key].domains[:-1]
        start = chunk[key].start_trim

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

    colourful_errors("INFO",
        "Parsing gene annotations.")

    file = variables.annotation_file_paths[0]

    if variables.annotation_type == "gff":
        basket, genes = gene_info_parser_gff(file)

    elif variables.annotation_type == "gb":
        basket, genes = gene_info_parser_genbank(file)
    
    if len(basket) == 0:
        colourful_errors("FATAL",
            "Watch out, no genomic features were loaded. The gb file is not being parsed correctly.")
        exit()

    contigs = {}
    for contig in variables.genome_seq:
        contigs[contig] = len(variables.genome_seq[contig])

    intergenic_regions = inter_gene_annotater(genes, contigs, variables.ir_size_cutoff)
    for ir in intergenic_regions:
        basket[ir] = Gene(gene=ir,
                            start=intergenic_regions[ir][0],
                            end=intergenic_regions[ir][1],
                            start_trim=int(intergenic_regions[ir][0]+(intergenic_regions[ir][2]*variables.subdomain_length[0])),
                            end_trim=int(intergenic_regions[ir][0]+(intergenic_regions[ir][2]*variables.subdomain_length[1])),
                            identity=ir,
                            contig=intergenic_regions[ir][3])
    return basket


def gene_basket_maker(basket,genes,gene_info):

    if gene_info["domain_size"] > variables.biggest_gene:
        variables.biggest_gene = gene_info["domain_size"]

    basket[gene_info["ID"]] = Gene(gene=gene_info["gene"], 
                                    start=gene_info["start"], 
                                    end=gene_info["end"], 
                                    orientation=gene_info["orientation"],
                                    identity=gene_info["ID"], 
                                    product=gene_info["product"], 
                                    contig=gene_info["contig"],
                                    start_trim=int(gene_info["start"]+(gene_info["domain_size"]*variables.subdomain_length[0])),
                                    end_trim=int(gene_info["start"]+(gene_info["domain_size"]*variables.subdomain_length[1])))

    genes.append((gene_info["start"],gene_info["end"],gene_info["orientation"],gene_info["gene"],gene_info["contig"]))
    return basket,genes

def gene_info_parser_genbank(file):

    basket = {}
    genes = []
    for rec in SeqIO.parse(file, "gb"):
        for feature in rec.features:
            gene_info = gb_parser(feature)
            if gene_info["ID"] != None:
                gene_info["contig"] = rec.id
                basket, genes = gene_basket_maker(basket,genes,gene_info)
    return basket, genes


def gene_info_parser_gff(file):

    genes = []
    basket = {}
    with open(file) as current:
        for line in current:
            line = line.split('\t')
            gene_info = gff_parser(line)
            if gene_info is None:
                break
            if gene_info:
                basket, genes = gene_basket_maker(basket,genes,gene_info)
    return basket, genes


def insertions_parser(startup=True):

    insertions_df = pd.read_csv(variables.insertion_file_path)
    insertions_df["unique"] = insertions_df["#Contig"] + \
        insertions_df["position"].astype(str)+insertions_df["Orientation"]
    insertions_df.drop_duplicates(subset=['unique'], inplace=True)
    if startup:
        colourful_errors("INFO",
            f"Total Insertions in library: {len(insertions_df)}")

    insertions_df.drop(
        ['Read Counts', 'Average mapQ across reads', "unique"], axis=1, inplace=True)

    insertions_df = insertions_df[insertions_df["Relative Position in Gene (0-1)"].astype(
        float) <= variables.subdomain_length[1]]
    insertions_df = insertions_df[insertions_df["Relative Position in Gene (0-1)"].astype(
        float) >= variables.subdomain_length[0]]

    variables.total_insertions = len(insertions_df)

    orient_pos = len([n for n in insertions_df["Orientation"] if n == "+"])
    variables.positive_strand_tn_ratio = orient_pos / variables.total_insertions

    for (contig, insertion, orientation, border) in zip(insertions_df["#Contig"],
                                                        insertions_df["position"],
                                                        insertions_df["Orientation"],
                                                        insertions_df["Transposon chromosome Border Sequence"]):

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
        borders = contig_df['Transposon chromosome Border Sequence'].tolist()
        
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
    
def insertion_annotater(chunk, variables):
    ''' The function insertion_annotater takes a dictionary chunk of gene 
    objects and a set of variables and annotates each gene object with 
    various insertion-related information, such as the number of insertions 
    and their orientation within each subdomain of the gene, the GC content of 
    each subdomain, and the transposon motif content of each subdomain. T
    he function returns the annotated chunk. '''

    for key in chunk:

        chunk[key].subdomain_insert_orient_neg, chunk[key].subdomain_insert_orient_plus = {}, {}  # clear past values
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
            domain_start = subdomain-chunk[key].start_trim
            domain_end = chunk[key].domains[i+1]-chunk[key].start_trim
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

    # divides the dictionary keys into smaller blocks that can be efficiently be multiprocessed
    divider = len(basket)//variables.cpus
    return_list = [dict() for i in range(variables.cpus)]
    i, list_iter = 0, 0
    for k in basket:
        if i < variables.cpus:
            return_list[i][k] = basket[k]

        if i == variables.cpus:  # odd number split will be distributed equally
            list_iter += 1
            return_list[list_iter][k] = basket[k]

        elif len(return_list[i]) >= divider:
            i += 1

    result_objs = []
    pool = multiprocessing.Pool(processes = variables.cpus)

    for chunk in return_list:
        result = pool.apply_async(
            insertion_annotater, args=((chunk, variables)))
        result_objs.append(result)

    pool.close()
    pool.join()
    result = [result.get() for result in result_objs]

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
    result_objs, pvalue_listing, euclidean_points = [], [], []

    pvaluing_array, names, pvalues_list = class_to_numba(basket)

    pool = multiprocessing.Pool(processes = variables.cpus)
    for p in pvalue:
        result = pool.apply_async(pvalue_iteration, args=((names, 
                                                           pvalues_list, 
                                                           pvaluing_array, 
                                                           p, 
                                                           pvalue_listing,
                                                           euclidean_points, 
                                                           variables)))
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

    for i, (_, entry) in enumerate(zip(names, pvaluing_array)):
        _, pvaluing_array = pvaluing_jit(
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
                              [int(optimal_basket[gene].significant[domain].dom_len)] +
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

    output_writer(variables.directory, f"{variables.output_name}_verbose", genes_list)

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

    stats_file = f"""\n    ####\nESSENTIALITY INFO\n    ####\n\nTransposon insertion bias (+strand, per dinucleotide, in %):\n{dic}\np-value cutoff: {pvalue}\nfdr corrected p-value cutoff: {fdr}\nNumber of features with at least one domain that is non-essential: {len(non_essentials_list)}\nNumber of whole features that are non-essential: {full_non_e_genes}\nNumber of features with at least one domain too small for assaying: {len(non_assayed_list)}\nNumber of whole features too small for assaying: {full_na_genes}\nNumber of features with at least one domain that is essential: {len(essentials)}\nNumber of whole features that are essential: {len(significant_genes_list)}\n"""
    text_out = variables.directory + "/Essentiality_stats.log"
    with open(text_out, "w+") as text_file:
        text_file.write(stats_file)

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
        f"{variables.directory}/ROC_curves_{variables.strain}.png", dpi=300)
    plt.close()
    gene_essentiality_compressor(optimal_basket)


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
                colourful_errors("INFO",
                    f"Loaded contig {contig}")
            variables.genome_length += len(variables.genome_seq[contig])

    elif variables.annotation_type == "gb":
        for rec in SeqIO.parse(variables.annotation_file_paths[0], "gb"):
            variables.genome_length = len(rec.seq)
            variables.annotation_contig = rec.id
            variables.borders_contig[variables.annotation_contig] = {}
            variables.orientation_contig[variables.annotation_contig] = {}
            variables.insertions_contig[variables.annotation_contig] = {}
            if startup:
                colourful_errors("INFO",
                    f"Loaded contig {variables.annotation_contig}")
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
            colourful_errors("INFO",
                f"Current domain division iteration: {i+1} out {iteration}")
            iterator_store = iterating(
                i, iterator_store, euclidean_distances, basket, current_gap)
            i += 1
            
        # sort by best eucledian distance
        sorted_optimal = sorted(iterator_store, key=lambda e: e[2])
        variables.best_domain_size = sorted_optimal[0][0]
        
    except Exception:
        colourful_errors("WARNING",
            "Low saturating library detected.")
        variables.best_domain_size = sorted_optimal[0]
    
    colourful_errors("INFO",
        f"Optimal domain division size of {variables.best_domain_size}bp")

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
    #output_writer(variables.directory, "Best_ROC_points{}".format(variables.output_name), iterator_store[best_index][4][0])

    return best_index, iterator_store[best_index][1], iterator_store[best_index][4]

def gene_essentiality_compressor(basket):
    
    essentials = pd.read_csv(f"{variables.directory}/{variables.output_name}_verbose.csv",
                             low_memory=False)

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
    main(sys.argv[0])
    multiprocessing.set_start_method("spawn")
