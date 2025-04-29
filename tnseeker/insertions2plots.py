import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns
from matplotlib.lines import Line2D
import sys
from Bio import SeqIO
from tnseeker.extras.helper_functions import variables_parser,adjust_spines,gb_parser,gff_parser

""" The script visualizes the distribution of insertions 
    in a genomic dataset. It processes the input data, processes genomic annotations, 
    and generates a plot with the distribution of insertions along the genome.                        
"""
    

def main(argv):
    global variables
    variables = variables_parser(argv)
    variables["barcode"] = variables["barcode"] == "True"
    gene_l,genome_seq = get_gene_len()
    plotter(genome_seq)
    reads_per_gene(gene_l)

def get_gene_len():
    
    gene_l = {}
    genome_seq={}
    if variables["annotation_type"] == "gb":
        for rec in SeqIO.parse(variables['annotation_file'], "gb"):
            genome_seq[rec.id] = len(rec.seq)
            for feature in rec.features:
                if feature.type != 'source':
                    gene_info = gb_parser(feature)
                    if gene_info:
                        if gene_info['gene'] not in gene_l:
                            gene_l[gene_info["gene"]] = gene_info["end"] - gene_info["start"]
                    
    else:
        with open(variables["annotation_file"]) as current:
            for line in current:
                line = line.split('\t')
                gene_info = gff_parser(line)
                if gene_info:
                    gene_l[gene_info["gene"]] = gene_info["end"] - gene_info["start"]

                if "sequence-region" in line[0]:
                    line = line[0].split(" ")
                    genome_seq[line[-3]] = int(line[-1][:-1])

    return gene_l,genome_seq

def reads_per_gene(gene_l):

    df = pd.read_csv(f"{variables['directory']}/all_insertions_{variables['strain']}.csv")
    dict_df = pd.DataFrame({'Gene Name': list(gene_l.keys()), 'lenght': list(gene_l.values())})
    gene_counts = df.groupby('Gene Name').size().reset_index(name='insertions')
    merged_df = pd.merge(gene_counts, dict_df, on='Gene Name', how='left')
    
    if variables["barcode"]:
        barcodes_per_gene(dict_df)
    
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
    #ax.set_xlim(0, 1)
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
    
def plotter(genome_seq):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    from matplotlib.lines import Line2D

    entry = []
    with open(f'{variables["directory"]}/all_insertions_{variables["strain"]}.csv') as current:
        for line in current:
            line = line.strip().split(",")
            if '#' not in line[0]:
                orientation = 0.82
                if line[2] == '+':
                    orientation = 0.25
                entry.append([line[0], int(line[1]), orientation, float(line[4])])

    entry = sorted(entry, key=lambda x: x[0])

    df = pd.DataFrame({
        'contig': [n[0] for n in entry],
        'position': [n[1] for n in entry],
        'orientation': [n[2] for n in entry],
        'reads': [n[3] for n in entry]
    })

    cmap = matplotlib.colormaps['jet']

    fig, ax = plt.subplots()  
    ax.remove()
    grid = plt.GridSpec(1, 5, wspace=0.4, hspace=0.3)
    ax=plt.subplot(grid[0, :4])

    incremental_position = 0
    genome_seq = {k: v for k, v in sorted(genome_seq.items(), reverse=True, key=lambda item: item[1])}

    max_y_value = 0
    bin_size = 25_000
    window = 5 

    for i, (contig, contig_length) in enumerate(genome_seq.items()):
        if contig in set(df['contig']):
            df2 = df[df['contig'] == contig].copy()

            df2['position'] += incremental_position
            df2['bin'] = (df2['position'] // bin_size) * bin_size

            grouped = df2.groupby('bin')
            bin_summary = grouped.agg(
                total_reads=('reads', 'sum'),
                total_insertions=('reads', 'count'),
                pos_reads=('reads', lambda x: x[df2.loc[x.index, 'orientation'] == 0.25].sum()),
                pos_insertions=('orientation', lambda x: (x == 0.25).sum()),
                neg_reads=('reads', lambda x: x[df2.loc[x.index, 'orientation'] == 0.82].sum()),
                neg_insertions=('orientation', lambda x: (x == 0.82).sum())
            ).reset_index()

            pos_insertions_avg = bin_summary['pos_insertions'].rolling(window=window, min_periods=1).mean()
            neg_insertions_avg = bin_summary['neg_insertions'].rolling(window=window, min_periods=1).mean()
            pos_reads_avg = bin_summary['pos_reads'].rolling(window=window, min_periods=1).mean()
            neg_reads_avg = bin_summary['neg_reads'].rolling(window=window, min_periods=1).mean()

            # Plot trailing averages
            ax.plot(bin_summary['bin'], pos_insertions_avg, linewidth=1.5, color=cmap(0.17), label='+ strand insertions (avg)')
            ax.plot(bin_summary['bin'], neg_insertions_avg, linewidth=1.5, color=cmap(0.91), label='- strand insertions (avg)')
            ax.plot(bin_summary['bin'], pos_reads_avg, linewidth=1.5, color=cmap(0.25), label='+ strand reads (avg)')
            ax.plot(bin_summary['bin'], neg_reads_avg, linewidth=1.5, color=cmap(0.82), label='- strand reads (avg)')
            
            ax.scatter(df2['position'], df2['reads'], c=cmap(df2['orientation']), s=0.3, alpha=0.4)

            max_y_value = max(
                max_y_value,
                df2['reads'].max(),
                bin_summary[['pos_reads', 'neg_reads']].max().max(),
                bin_summary[['pos_insertions', 'neg_insertions']].max().max()
            )

        # Background shading for contigs
        current_start = incremental_position
        current_end = incremental_position + contig_length
        if i % 2 == 0:
            ax.axvspan(current_start, current_end, facecolor='lightgrey', alpha=0.3)
        else:
            ax.axvspan(current_start, current_end, facecolor='white', alpha=0.1)

        incremental_position += contig_length

    # Axis styling
    ax.set_yscale('log')
    ax.set_ylabel("Amount of reads or insertions")
    ax.set_xlabel("Cumulative genome position (bp)")
    ax.set_title(variables["strain"])

    adjust_spines(ax, ['left', 'bottom'], (0, ax.get_xlim()[1]), (1,max_y_value))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    legend_elements = [
        Line2D([0], [0], color=cmap(0.25), label='+ strand reads'),
        Line2D([0], [0], color=cmap(0.17), label='+ strand insertions'),
        Line2D([0], [0], color=cmap(0.82), label='- strand reads'),
        Line2D([0], [0], color=cmap(0.91), label='- strand insertions')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=8)

    ax.set_ylim(1, max_y_value * 1.1)
    ax.set_xlim(0, incremental_position)

    # Side KDE plot
    ax2=plt.subplot(grid[0, 4])
    ax2.set_ylim(0, np.log(max_y_value))

    plt.yticks([])

    sns.kdeplot(
               y=np.log(df['reads']),
               fill=True, common_norm=False,
               alpha=.5, linewidth=0)

    adjust_spines(ax2, ['left', 'bottom'], (0, ax2.get_xlim()[1]), (0, np.log(max_y_value)))
    

    plt.subplots_adjust(bottom=0.15)  # Ensure legend/axis labels don't get clipped
    plt.savefig(f'{variables["directory"]}/reads_insertion_genome_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main(sys.argv[0])
