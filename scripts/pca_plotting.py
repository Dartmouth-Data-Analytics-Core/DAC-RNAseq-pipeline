#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
import matplotlib.patches as mpatches

import sys

# Argument ordering (for sys)
# sys[0] = script name (not an argument)
# sys[1] = input tsv feature counts file
# sys[2] = the number of genes with largest variance to create PCA vectors from
# sys[3] = Path name for output plots (folder name, not including final "/")
# sys[4] = (optional, default = 2) Number of Highest PCA genes to compare
#          and pairwise plot points on (must include if including sys[5])
# sys[5] = (optional) path to tsv file containing hex values for sample groups


#df = pd.read_csv("test_gene_counts.tsv", sep='\t', index_col=0)
df = pd.read_csv(sys.argv[1], sep='\t', index_col=0)
#df = pd.read_csv("../featurecounts.readcounts.tsv", sep='\t', index_col=0)
#Remove featurecounts info columns, if necessary. Will be necessary for pipeline

try:
    df = df.drop(["Gene Name","Chr","Start","End","Strand","Length"], axis=1)
except KeyError:
    pass
try:
    df = df.drop(["Chr","Start","End","Strand","Length"], axis=1)
except KeyError:
    pass



sample_names = list(df.columns)

### Transform data
#log2 scale dataframe

def log2_normalize(df):
    x = None
    try:
        x = np.log2(df)
    except:
        x = np.log2((df + 1))
    return x
def ln_normalize(df):
    x = None
    try:
        x = np.ln(df)
    except:
        x = np.ln((df + 1))
    return x

df = log2_normalize(df)


############# Load Metadata (For Colors and Type Labels)

def extract_metadata(file_location):
    groups_df = pd.read_csv(file_location, sep='\t', index_col=0)
    groups_df = groups_df.applymap(str)
    # groups_df = groups_df.T[sample_names].T
    assert(len(sample_names) == len(groups_df.index))
    types = groups_df['group'].unique()
    types.sort()
    # If "Control" or "control" is a label, bring to top
    lft_lower = types.searchsorted("control", side = 'left')
    rght_lower = types.searchsorted("control", side = 'right')
    lft_Upper = types.searchsorted("Control", side = 'left')
    rght_Upper = types.searchsorted("Control", side = 'right')
    types_list = types.tolist()

    if lft_lower != rght_lower:
        lft_Upper += 1
        rght_Upper += 1
        types_list.pop(lft_lower)
        types_list.insert(0, "control")


    if lft_Upper != rght_Upper:
        types_list.pop(lft_Upper)
        types_list.insert(0, "Control")
    groups = [list(np.where(groups_df['group'] == types_list[i])[0]) for i in range(len(types_list))]
    # while len(groups[-1]) == 0:
    #     groups.pop()
    return groups_df, types_list, groups



groups_df = pd.DataFrame()
groups = None

try:
    groups_df, types_list, groups = extract_metadata("test_metadata.tsv")

except:
    try: 
        groups_df, types_list, groups = extract_metadata("sample_fastq_list.txt")

    except:
        try:
            groups_df, types_list, groups = extract_metadata("sample_fastq_paired_list.txt")

        except:
            groups = [[i for i in range(len(sample_names))]]
            # if metadata doesn't exist, filler dataframe only
            # for notification purposes that metadata doesn't exist.
            # can be improved                                                                                        
            groups_df = pd.DataFrame()
    

# group_1 = list(np.where(groups_df['group'] == "1")[0])
# group_2 = list(np.where(groups_df['group'] == "2")[0])
# group_3 = list(np.where(groups_df['group'] == "3")[0])
# group_4 = list(np.where(groups_df['group'] == "4")[0])
# group_5 = list(np.where(groups_df['group'] == "5")[0])
# group_6 = list(np.where(groups_df['group'] == "6")[0])

#### Set up coloring
color_list = ['#4DBBD5', '#E64B35', '#00A087', '#8473E2', '#EAA61D', '#D069D4']
color_list_temp = ['#808080' for i in range(max(len(groups), 0))]
color_list += color_list_temp
color_list = color_list[:len(groups)]

# if there is a file specifying colors, load the colors
try:
    df_temp = None
    if len(sys.argv) > 5:
        df_temp = pd.read_csv(sys.argv[5], sep='\t', index_col=0)
    else:
        df_temp = pd.read_csv("sample_colors_hex.tsv", sep='\t', index_col=0)
    for i in range(len(groups)):
        color_list[i] = df_temp['color'][i]
except FileNotFoundError:
    pass


group_colors = []


# If no metadata was found, color all points the same. Otherwise, color by sample type
if groups_df.shape[0] != 0:
    for i in range(len(groups_df['group'])):
        group_colors.append(color_list[types_list.index(groups_df['group'][i]) - 1])
else:
    group_colors = [color_list[0] for i in range(len(sample_names))]
  
          
    # group_colors.append(color_list[])
    #     group_colors.append('#E64B35')
    # if i in group_3:
    #     group_colors.append('#00A087')
    # if i in group_4:
    #     group_colors.append('#8473E2')     
    # if i in group_5:
    #     group_colors.append('#EAA61D')
    # if i in group_6:
    #     group_colors.append('#D069D4')


################################################################

#number_of_top_genes = 500
number_of_top_genes = int(sys.argv[2])

def gene_variance_filter(df, number_of_top_var_genes = 500):
### Select high variance genes
    #Using Variance
    gene_var = df.var(axis=1) 
    #Using Coefficient of Variation
    # gene_var = df.var(axis=1)**2 / df.mean(axis=1)
    sorted_var = gene_var.sort_values(ascending=False)
    top_var = sorted_var.head(number_of_top_var_genes)
    genes_to_consider = top_var.index.values.tolist()
    df = df.loc[genes_to_consider]
    return df

df = gene_variance_filter(df, number_of_top_genes)



#### Heatmap Plotting ####    
                
def plot_heatmap(path_name = "."):
    cg = sns.clustermap(df, col_colors=group_colors,  cmap="vlag", vmin=-3, vmax=3, z_score=0, center=0, xticklabels=True, yticklabels=False, dendrogram_ratio=.085, cbar_pos=(0.005, .878, .04,.115), figsize=(12, 18))
    # cg = sns.clustermap(df,  cmap="vlag", xticklabels=True, yticklabels=False, dendrogram_ratio=.085, cbar_pos=(0.005, .878, .04,.115), figsize=(12, 18))
    
    cg.ax_col_dendrogram.set_visible(True)
    ax = cg.ax_heatmap
    
    legpatch = [mpatches.Patch(color=color_list[i], label=types_list[i]) for i in range(len(color_list))]
    # legpatch2 = mpatches.Patch(color='#E64B35', label="2")
    # legpatch3 = mpatches.Patch(color='#00A087', label="3")                           
    # legpatch4 = mpatches.Patch(color='#8473E2', label="4")
    
    cg.ax_col_dendrogram.legend(handles=legpatch, bbox_to_anchor=(.98,.99),
                bbox_transform=plt.gcf().transFigure, fontsize=12, title="Group", title_fontsize="x-large")
    
    cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_ymajorticklabels(), fontsize = 14)
    cg.savefig(path_name + "/Heatmap_scaled_" + str(number_of_top_genes)+ "_features.png")

plot_heatmap(sys.argv[3])


#### PCA Plotting ####
from sklearn.preprocessing import quantile_transform

def PCA_plot(path_name = ".", num_PCA_vec_plots = 2, quantile_normalize = False, include_sample_names = True, plot_length = 22/2.54, plot_width = 16/2.54):
    #optional normalization before PCA

    if quantile_normalize:
        df_new = quantile_transform(df, axis=1)
    else:
        df_new = df
    #Compute PCA
    pca = PCA(n_components=num_PCA_vec_plots, svd_solver='full')
    principalComponents = pca.fit_transform(df_new.T)
    pcadf = pd.DataFrame(data = principalComponents)
    
    for more_sig in range(num_PCA_vec_plots): # pca vector with more variance, so smaller number
        for less_sig in range(more_sig + 1, num_PCA_vec_plots): # pca vector with less variance, so larger number
            
            plt.figure(figsize=(plot_length, plot_width))
            
            #Set up Legend
            
            legpatch = [mpatches.Patch(color=color_list[i], label=types_list[i]) for i in range(len(groups))]
            plt.legend(handles=legpatch, bbox_to_anchor=(.87,.92),
                        bbox_transform=plt.gcf().transFigure, fontsize=14, title="Group", title_fontsize="x-large")
            
            #Draw scatterplot points
            for i in range(len(groups)):
                plt.scatter(pcadf[more_sig].iloc[groups[i]],(pcadf[less_sig]*1).iloc[groups[i]], s=240, alpha=.85, edgecolors='black', linewidth=1, marker='o', c=color_list[i])
            
            
            #Add sample names
            if include_sample_names:
                for i, txt in enumerate(sample_names):
                    plt.text(pcadf[more_sig].iloc[i]+0, pcadf[less_sig].iloc[i]+0, txt, size=16)
            
            #Label X and Y axes
            plt.xlabel("PC" + str(more_sig + 1) + ": "+ str(int(np.round(pca.explained_variance_ratio_[more_sig], decimals=2)*100))  +"% variance explained", size=16)
            plt.ylabel("PC" + str(less_sig + 1) + ": "+ str(int(np.round(pca.explained_variance_ratio_[less_sig], decimals=2)*100))  +"% variance explained", size=16)
            
            #Resize axis tick labels
            ax = plt.gca()
            ax.set_xticklabels(ax.get_xticks(), fontsize=14)
            ax.set_yticklabels(ax.get_yticks(), fontsize=14)
            
            #Remove the top and right frame from the plot.
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            #Tighten Layout
            plt.tight_layout()
            plt.savefig(path_name + '/PCA_' + str(more_sig + 1) + '_vs_' + str(less_sig + 1) + '.png')

if len(sys.argv) > 4:
    PCA_plot(sys.argv[3], int(sys.argv[4]))
else:
    PCA_plot(sys.argv[3], 2)


###  PCA Variance Ratio Bar Plot  ####

def PCA_Variance_Plot(path_name = ".", num_PCA_vect = 6, quantile_normalize = False, plot_length = 22/2.54, plot_width = 16/2.54):
    
    if quantile_normalize:
        df_new = quantile_transform(df, axis=1)
    else:
        df_new = df
    
    pca = PCA(n_components=num_PCA_vect, svd_solver='full')
    principalComponents = pca.fit_transform(df_new.T)
    pcadf = pd.DataFrame(data = principalComponents)
    
    plt.figure(figsize=(plot_length, plot_width))
    plt.title("PCA Variance Ratios")
    # Define the bar parameter lists
    plt.xlim(0.5, num_PCA_vect + 0.5)
    bar_x_coords = [(i+1) for i in range(num_PCA_vect)]
    bar_y_tops = [pca.explained_variance_ratio_[i] for i in range(num_PCA_vect)]
    labels = ["PC"+str(i+1) for i in range(num_PCA_vect)]
    # bar_y_bottoms = [0 for i in range(num_PCA_vect)]
    
    # Create the Plot
    plt.bar(bar_x_coords, bar_y_tops, tick_label=labels)
    
    #Tighten Layout
    plt.tight_layout()
    plt.savefig(path_name + '/PCA_Variance_Bar_Plot.png')
    
PCA_Variance_Plot(sys.argv[3])