import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns
from matplotlib.lines import Line2D
import sys
from Bio import SeqIO
from tnseeker.extras.helper_functions import variables_parser,adjust_spines

""" The script visualizes the distribution of insertions 
    in a genomic dataset. It processes the input data, processes genomic annotations, 
    and generates a plot with the distribution of insertions along the genome. 
    Here's a step-by-step breakdown of the code:
                                
"""
    

def main(argv):
    global variables
    variables = variables_parser(argv)
    variables["barcode"] = variables["barcode"] == "True"
    plotter()
    reads_per_gene()

def get_gene_len():
    gene_l = {}
    
    if variables["annotation_type"] == "gb":
        for rec in SeqIO.parse(variables["annotation_file"], "gb"):
            for feature in rec.features:
                if (feature.type != 'source') & (feature.type != 'misc_feature'):
                    start = feature.location.start
                    end = feature.location.end
                    domain_size = end - start
                    
                    try:
                        identity = feature.qualifiers['locus_tag'][0]
                        if "gene" in feature.qualifiers:
                            identity = feature.qualifiers['gene'][0]
                        gene_l[identity] = domain_size
                        
                    except KeyError:
                        identity = None
                    
    else:
        with open(variables["annotation_file"]) as current:
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
                    try:
                        if "gene" in feature:
                            gene=feature["gene"]
                        if "Name" in feature:
                            gene=feature["Name"]
                        else:
                            gene=feature["ID"]
                    
                    except KeyError:
                        identity = None
                    
                    gene_l[gene] = domain_size
                    
    return gene_l

def reads_per_gene():

    gene_l = get_gene_len()
    
    df = pd.read_csv(f"{variables['directory']}/all_insertions_{variables['strain']}.csv")
    dict_df = pd.DataFrame({'Gene Name': list(gene_l.keys()), 'lenght': list(gene_l.values())})
    gene_counts = df.groupby('Gene Name').size().reset_index(name='insertions')
    merged_df = pd.merge(gene_counts, dict_df, on='Gene Name', how='left')
    
    if variables["barcode"]:
        barcodes_per_gene(dict_df,variables["directory"],variables["strain"])
    
    merged_df["insertions/gene_len"] = merged_df['insertions'] / merged_df['lenght']
    merged_df = merged_df.dropna(subset=['insertions/gene_len'])
    
    fig, ax = plt.subplots()  
    sns.histplot(
        data=merged_df,
        x='insertions/gene_len',
        kde=False
    )
    
    plt.ylabel('Number of genes')
    plt.xlabel('Tn insertions per gene per length')
    plt.title('Histogram of transposon insertions per gene per length')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    adjust_spines(ax, ['left', 'bottom'], (0, ax.get_xlim()[1]), (1,ax.get_ylim()[1]))
    plt.savefig(f'{variables["directory"]}/Histogram of transposon insertions per gene per length.png', dpi=300,bbox_inches='tight')
    
    plt.show()

def barcodes_per_gene(dict_df):
    df = pd.read_csv(f'{variables["directory"]}/annotated_barcodes_{variables["strain"]}.csv')
    gene_counts = df.groupby('Gene Name').size().reset_index(name='#Barcode')
    
    merged_df = pd.merge(gene_counts, dict_df, on='Gene Name', how='left')
    merged_df["barcodes/gene_len"] = merged_df['#Barcode'] / merged_df['lenght']
    merged_df = merged_df.dropna(subset=['barcodes/gene_len'])
    
    fig, ax = plt.subplots()  
    sns.histplot(
        data=merged_df,
        x='barcodes/gene_len',
        kde=False
    )
    
    plt.ylabel('Number of genes')
    plt.xlabel('barcodes per gene per length')
    plt.title('Histogram of barcodes per gene per length')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    adjust_spines(ax, ['left', 'bottom'], (0, ax.get_xlim()[1]), (1,ax.get_ylim()[1]))
    plt.savefig(f'{variables["directory"]}/Histogram of barcodes per gene per length.png', dpi=300,bbox_inches='tight')
    
    plt.show()
    
def plotter():
    entry=[]
    with open(f'{variables["directory"]}/all_insertions_{variables["strain"]}.csv') as current:
        for line in current:
            line=line[:-1].split(",")
            if '#' not in line[0]:
                c=0.82
                if line[2] == '+':
                    c=0.25
                entry.append([line[0]]+[int(line[1])]+[c]+[float(line[4])])
    
    entry = sorted(entry,key=lambda x: x[0])
    
    df = pd.DataFrame({'contig':[n[0] for n in entry],
                      'position':[n[1] for n in entry],
                      'orientation':[n[2] for n in entry],
                      'reads':[n[3] for n in entry]})
    
    contig_max = {}
    for entry in df['contig']:
        if entry not in contig_max:
            contig_max[entry] = max(df[df['contig']==entry]['position'])

    genome_seq={}
    if variables["annotation_type"] != 'gb':
        contig=None
        with open(variables["genome_file"]) as current: #NT12004_22
            for line in current:
                if '>' not in line:
                    line = line[:-1]
                    genome_seq[contig] += len(line)
                else:
                    contig = line[1:-1]
                    genome_seq[contig] = 0
    else:
        for rec in SeqIO.parse(variables["genome_file"], "gb"):
            genome_seq[rec.id] = len(rec.seq)
          
    cmap = matplotlib.colormaps['jet'] #gnuplot
    colour = ['.9','.6']*len(genome_seq)
    
    fig, ax = plt.subplots()  
    ax.remove()
    grid = plt.GridSpec(1, 5, wspace=0.4, hspace=0.3)
    ax=plt.subplot(grid[0, :4])
    max_y_value_cum_pos,max_y_value_cum_neg,max_y_value_reads,max_y_value_pos,max_y_value_neg = 0,0,0,0,0
    incremental_position=0
    genome_seq={k: v for k, v in sorted(genome_seq.items(), reverse=True,key=lambda item: item[1])}
    
    for i,contig in enumerate(genome_seq):

        if contig in df['contig'].values:
            df2 = df[df['contig']==contig]
            if i>0:
                df2.loc[:, 'position'] += incremental_position #sum of previous positions into the new contig as if it were concatenated
            
            plt.scatter(x=df2['position'], 
                        y=df2['reads'],
                        c=cmap(df2['orientation']),
                        s=.3,
                        alpha=.4)
            max_y_value_reads = max(max_y_value_reads, df2['reads'].max())
            
            ## binning the reads per genome position
            k,inserts,inserts_pos,inserts_neg=100000+incremental_position,0,0,0
            inserts_pos_cum,inserts_neg_cum=0,0
            binned,binned_pos,binned_neg,binned_pos_cum,binned_neg_cum = {},{},{},{},{}
            df2 = df2.sort_values('position')
            for j,entry in enumerate(df2['position']):
                inserts += 1
                if df2['orientation'].iloc[j] == 0.25:
                    inserts_pos+=df2['reads'].iloc[j]
                    inserts_pos_cum+=1
                else:
                    inserts_neg+=df2['reads'].iloc[j]
                    inserts_neg_cum+=1
                if entry>k:
                    binned[k] = inserts
                    binned_pos_cum[k] = inserts_pos_cum
                    binned_neg_cum[k] = inserts_neg_cum
                    binned_pos[k] = inserts_pos
                    binned_neg[k] = inserts_neg
                    inserts,inserts_pos,inserts_neg,inserts_pos_cum,inserts_neg_cum=0,0,0,0,0
                    k+=100000
            
            df3=pd.DataFrame.from_dict([binned_pos_cum]).transpose()
            plt.plot(df3.index,
                     df3[0],
                     linewidth=1.5,
                     color=cmap(0.17))
            max_y_value_cum_pos = max(max_y_value_cum_pos, df3[0].max())
            
            df3=pd.DataFrame.from_dict([binned_neg_cum]).transpose()
            plt.plot(df3.index,
                     df3[0],
                     linewidth=1.5,
                     color=cmap(0.91))
            max_y_value_cum_neg = max(max_y_value_cum_neg, df3[0].max())
            
            pos=pd.DataFrame.from_dict([binned_pos]).transpose()
            plt.plot(pos.index,
                     pos[0],
                     linewidth=1.5,
                     color=cmap(0.25)) #0.25
            max_y_value_pos = max(max_y_value_pos, pos[0].max())
            
            neg=pd.DataFrame.from_dict([binned_neg]).transpose()
            plt.plot(neg.index,
                     neg[0],
                     linewidth=1.5,
                     color=cmap(0.82)) #0.82
            max_y_value_neg = max(max_y_value_neg, neg[0].max())

        incremental_position+=genome_seq[contig]
        
        #ax.set_ylim(0.95, 1.2)
        ax.set_yscale('log')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        
        plt.ylabel("Amount of reads or insertions")
        plt.xlabel("Cumulative genome position (bp)")
        plt.title(variables["strain"])
        
        max_y_value = max([max_y_value_cum_pos,max_y_value_cum_neg,
                           max_y_value_reads,
                           max_y_value_pos,max_y_value_neg])
        max_y_value *= 1.1
        
        adjust_spines(ax, ['left', 'bottom'], (0, ax.get_xlim()[1]), (1,max_y_value))
    
    for i,contig in enumerate(genome_seq):
        incremental_position=0
        plt.fill_between([incremental_position,genome_seq[contig]+incremental_position],0,max_y_value,alpha=0.2,color=colour[i])
        incremental_position+=genome_seq[contig]
    
    legend_elements = [Line2D([0], [0], color=cmap(0.25), label='+ strand reads'),
                       Line2D([0], [0], color=cmap(0.17), label='+ strand insertions'),
                       Line2D([0], [0], color=cmap(0.82), label='- strand reads'),
                       Line2D([0], [0], color=cmap(0.91), label='- strand insertions')]
                               
    ax.legend(handles=legend_elements, loc='upper center',bbox_to_anchor=(0.5, -0.15),ncol=2,prop={'size': 8})
    ax.set_ylim(1, max_y_value)
    
    ax2=plt.subplot(grid[0, 4])
    ax2.set_ylim(0, np.log(max_y_value))

    plt.yticks([])

    sns.kdeplot(
               y=np.log(df['reads']),
               fill=True, common_norm=False,
               alpha=.5, linewidth=0)

    adjust_spines(ax2, ['left', 'bottom'], (0, ax2.get_xlim()[1]), (0, np.log(max_y_value)))
    
    plt.savefig(f'{variables["directory"]}/reads.png', dpi=300,bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main(sys.argv[0])
