import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns
from matplotlib.lines import Line2D
import sys
from Bio import SeqIO

""" The code is a Python script that visualizes the distribution of insertions 
    in a genomic dataset. It processes the input data, processes genomic annotations, 
    and generates a plot with the distribution of insertions along the genome. 
    Here's a step-by-step breakdown of the code:
    
    Define the plotter function that takes the directory, annotation file, 
    annotation type, and strain as input arguments and does the following:
        
    1. Read the input file containing the insertions data and store the relevant 
    information in a list called 'entry'.
    
    2. Sort the 'entry' list based on the contig and create a pandas DataFrame from it.
    
    3. Process the annotation file to obtain the length of each contig in the genome.
    
    4. Create a color map using the 'jet' color scheme from matplotlib.
    
    5. Initialize a figure and axis for plotting, along with setting up grid specifications.
    
    6. Iterate through the contigs in the genome and create a scatter plot showing 
    the distribution of insertions (with separate colors for the positive and negative strands).
    
    7. Perform binning of reads per genome position and plot them on the same axis.
    
    8. Customize the appearance of the axis, labels, and title of the plot.
    
    9. Add a legend to the plot.
    
    10. Create a second axis to display the distribution of reads as a histogram.
    
    11. Save the generated plot as an image file (reads.png) in the given directory.
    
    The script is executed with command-line arguments. It checks if there are 
    at least two command-line arguments and calls the main function with the 
    provided arguments."""
    

def main(argv):
    directory = argv[0]
    annotation = argv[1]
    anno_type = argv[2]
    strain = argv[3]

    plotter(directory,annotation,anno_type,strain)

def adjust_spines(ax, spines,x,y): #offset spines
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 2.5))  # outward by 10 points
            if (loc == 'left') or (loc == 'right'):
                spine.set_bounds(y)
            else:
                spine.set_bounds(x)
        else:
            spine.set_color('none')  # don't draw spine

def plotter(directory,annotation,anno_type,strain):
    entry=[]
    with open(f"{directory}/all_insertions_{strain}.csv") as current:
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
    if anno_type != 'gb':
        contig=None
        with open(annotation) as current: #NT12004_22
            for line in current:
                if '>' not in line:
                    line = line[:-1]
                    genome_seq[contig] += len(line)
                else:
                    contig = line[1:-1]
                    genome_seq[contig] = 0
    else:
        for rec in SeqIO.parse(annotation, "gb"):
            genome_seq[rec.id] = len(rec.seq)
          
    cmap = matplotlib.cm.get_cmap('jet') #gnuplot
    colour = ['.9','.6']*len(genome_seq)
    
    fig, ax = plt.subplots()  
    grid = plt.GridSpec(1, 5, wspace=0.4, hspace=0.3)
    ax=plt.subplot(grid[0, :4])
    ax.set_ylim(1, 1000000)
    incremental_position=0
    genome_seq={k: v for k, v in sorted(genome_seq.items(), reverse=True,key=lambda item: item[1])}
    
    for i,contig in enumerate(genome_seq):
        
        #plt.plot([incremental_position,genome_seq[contig]+incremental_position],[10,10],linewidth=5)
        plt.fill_between([incremental_position,genome_seq[contig]+incremental_position],0,1000000,alpha=0.2,color=colour[i])
        
        if contig in df['contig'].values:
            df2 = df[df['contig']==contig]
            if i>0:
                df2['position'] += incremental_position #sum of previous positions into the new contig as if it were concatenated
            
            plt.scatter(x=df2['position'], 
                        y=df2['reads'],
                        c=cmap(df2['orientation']),
                        s=.3)
            
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
            
            df3=pd.DataFrame.from_dict([binned_neg_cum]).transpose()
            plt.plot(df3.index,
                     df3[0],
                     linewidth=1.5,
                     color=cmap(0.91))
            
            pos=pd.DataFrame.from_dict([binned_pos]).transpose()
            plt.plot(pos.index,
                     pos[0],
                     linewidth=1.5,
                     color=cmap(0.25)) #0.25
            
            neg=pd.DataFrame.from_dict([binned_neg]).transpose()
            plt.plot(neg.index,
                     neg[0],
                     linewidth=1.5,
                     color=cmap(0.82)) #0.82

        incremental_position+=genome_seq[contig]
        
        #ax.set_ylim(0.95, 1.2)
        ax.set_yscale('log')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        
        plt.ylabel("Reads / insertions")
        plt.xlabel("Cumulative genome position (bp)")
        plt.title(f"{strain}")
        
        adjust_spines(ax, ['left', 'bottom'], (0, ax.get_xlim()[1]), (ax.get_ylim()))
    
    legend_elements = [Line2D([0], [0], color=cmap(0.25), label='+ strand reads'),
                       Line2D([0], [0], color=cmap(0.17), label='+ strand insertions'),
                       Line2D([0], [0], color=cmap(0.82), label='- strand reads'),
                       Line2D([0], [0], color=cmap(0.91), label='- strand insertions')]
                               
    ax.legend(handles=legend_elements, loc='upper center',bbox_to_anchor=(0.5, -0.15),ncol=2,prop={'size': 8})
    
    ax2=plt.subplot(grid[0, 4])
    ax2.set_ylim(0, np.log(1000000))

    plt.yticks([])

    sns.kdeplot(
               y=np.log(df['reads']),
               fill=True, common_norm=False,
               alpha=.5, linewidth=0)

    adjust_spines(ax2, ['left', 'bottom'], (0, ax2.get_xlim()[1]), (ax2.get_ylim()))
    
    plt.savefig(f"{directory}/reads.png", dpi=300,bbox_inches='tight')

if __name__ == "__main__":
    if len(sys.argv) > 2:
        argv = [sys.argv[1]]+[sys.argv[2]]+[sys.argv[3]]+[sys.argv[4]]
    main(argv)
