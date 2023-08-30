import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Make figure texts editable in Adobe Illustrator
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42

pd.set_option("display.max_columns", 50)

### Data reading helper functions ###
def read_expression_data(path):
    print(f"Reading expression data from {os.path.join(path, 'expression')}")
    counts = pd.read_csv(os.path.join(path, "expression", "exon_count_matrix.txt.gz"),
                    sep = '\t', index_col = 0)

    rpkms  = pd.read_csv(os.path.join(path, "expression" ,"exon_rpkm_matrix.txt.gz"),
                    sep = '\t', index_col = 0)

    print("Counts shape: ", counts.shape)
    print("RPKMs shape:  ", rpkms.shape)
    
    if counts.shape == rpkms.shape:
        print("The shape of the expression matrices match")
    else:
        print("WARNING: The shape of the expression matrices do not match")
        
    print("\n")
    return counts, rpkms

def read_meta(path):
    print(f"Reading meta data from {os.path.join(path, 'meta', 'meta.txt')}")
    meta = pd.read_csv(os.path.join(path, "meta", "meta.txt"), sep = '\t')
    meta.fillna('empty', inplace = True)
    print("\n")
    return meta

def read_conversion_rates(path, meta):
    print(f"Reading conversion rates from {os.path.join(path, 'stats', 'conversionrates.csv')}")
    conv_rates = pd.read_csv(os.path.join(path, "stats", 'conversionrates.csv'), index_col = 0)
    conv_rates = conv_rates[~conv_rates.conversionType.str.contains('n', case = False)]
    conv_rates = pd.merge(conv_rates, meta, left_on = 'cell', right_on = 'BC', how = 'inner')
    print("\n")
    return conv_rates

def read_pigs(path, CI_threshold = None):
    if CI_threshold:
        print(f"Reading pi_g values from {os.path.join(path, 'pigs' ,'pigs.txt')}\nApplying a C.I width threshold of {CI_threshold}")
        pigs = pd.read_csv(os.path.join(path, 'pigs' ,"pigs.txt"), sep = '\t', index_col = 0).query(f"conf_width <= {CI_threshold}")
        
    else:
        print(f"Reading pi_g values from {os.path.join(path, 'pigs' ,'pigs.txt')}\nNo C.I threshold supplied")
        pigs = pd.read_csv(os.path.join(path, 'pigs' ,"pigs.txt"), sep = '\t', index_col = 0)

    pig_matrix = pd.pivot_table(values = 'pi_g', index = 'gene_id', columns = 'BC', data = pigs).rename_axis(None, axis = 1).rename_axis(None, axis = 0)
    # Set NA-values to zero. They arise from not-detected genes in certain conditions.
    pig_matrix = pig_matrix.rename_axis(None, axis = 1).rename_axis(None, axis = 0)
    pig_matrix = pig_matrix.fillna(0)
    print("\n")    
    return pigs, pig_matrix

def read_conversion_probabilities(path):
    print(f"Reading conversion probabilities from {os.path.join(path, 'pcs', 'conv_probs.txt')}")
    conv_probs = pd.read_csv(os.path.join(path, "pcs", "conv_probs.txt"), sep = '\t')
    conv_probs.loc[:, 'snr'] = conv_probs.pc / conv_probs.pe
    print("\n")
    return conv_probs

def read_sequencing_stats(path, meta, pigs, conv_probs):
    print(f"Reading sequencing stats from {os.path.join(path, 'stats', 'readspercell.txt')}")
    stats_df = pd.read_csv(os.path.join(path, "stats", "readspercell.txt"), sep = '\t')
    stats_df = stats_df[stats_df.BC.isin(meta.BC.values)]
    
    # All features
    genes_with_pigs = pigs.groupby('BC')['gene_id'].nunique()
    
    # Genes with pi_g > 0
    genes_with_pigs_gtzero = pigs.query('pi_g > 0').groupby('BC').count().iloc[:, [0]].rename({'gene_id': 'num_features'}, axis = 1)
    
    stats_df = pd.merge(stats_df, genes_with_pigs, left_on = 'BC', right_index = True).rename({'gene_id': 'ngenes_with_pigs'}, axis = 1)
    stats_df = pd.merge(stats_df, conv_probs.loc[:, ['cell', 'snr']], left_on = 'BC', right_on = 'cell').drop('cell', axis = 1)
    stats_df = pd.merge(stats_df, genes_with_pigs_gtzero, left_on = 'BC', right_index = True).rename({"num_features": "ngenes_with_pigs_gtzero"}, axis = 1)
    
    print("\n")
    return stats_df

def save_conditions_to_keep(stats_df, outpath = None, readcount_threshold = 200_000, ngene_threshold = 7000, exon_threshold = 45, snr_threshold = 10, ngene_accessor = 'ngenes_with_pigs'):
    
    if not ngene_accessor in ["ngenes_with_pigs", "ngenes_with_pigs_gtzero"]:
        raise ValueError("Error: ngene_accessor needs to be one of [ngenes_with_pigs/ ngenes_with_pigs_gtzero]")
        
    conds_to_keep = stats_df.query(f"(Exon / total_reads)*100 > {exon_threshold} \
    and total_reads > {readcount_threshold} \
    and {ngene_accessor} > {ngene_threshold}\
    and snr > {snr_threshold}").BC.values.tolist()
    
    if outpath:
        if not os.path.exists(os.path.split(outpath)[0]):
            os.makedirs(os.path.split(outpath)[0])
        with open(os.path.join(outpath), 'w') as fh:
            print(f"Writing kept conditions to {outpath}")
            for bc in conds_to_keep:
                fh.write(f"{bc}\n")
        return conds_to_keep

    else:
        return conds_to_keep
    
def filter_conditions(counts, rpkms, pig_matrix, meta, stats_df, conds_to_keep):
    print("Filtering conditions...")
    
    # Filter expression matrices
    filtered_counts  = counts.loc[:, counts.columns.isin(conds_to_keep)]
    filtered_rpkms   = rpkms.loc[:, rpkms.columns.isin(conds_to_keep)]
    filtered_pig_matrix  = pig_matrix.loc[:, pig_matrix.columns.isin(conds_to_keep)]
    
    # Filter meta table and stats df
    filtered_meta      = meta[meta.BC.isin(conds_to_keep)]
    filtered_stats_df  = stats_df[stats_df.BC.isin(conds_to_keep)]
    
    filtered_counts = filtered_counts.sort_index(axis = 1).sort_index(axis = 0)
    filtered_rpkms  = filtered_rpkms.sort_index(axis = 1).sort_index(axis = 0)
    filtered_pig_matrix  = filtered_pig_matrix.sort_index(axis = 1).sort_index(axis = 0)
    
    # Use common genes
    common_genes = np.intersect1d(filtered_counts.index.values, filtered_pig_matrix.index.values)
    filtered_counts = filtered_counts.loc[common_genes, :]
    filtered_rpkms  = filtered_rpkms.loc[common_genes, :]
    filtered_pig_matrix  = filtered_pig_matrix.loc[common_genes, :]
    
    return filtered_counts, filtered_rpkms, filtered_pig_matrix, filtered_meta, filtered_stats_df


def separate_expression_matrices(counts, rpkms, pig_matrix, outfolder = None):
    print("Separating expression matrices into new and pre-existing RNA...")
    
    new_counts = counts * pig_matrix
    old_counts = counts * (1-pig_matrix)

    new_rpkms = rpkms * pig_matrix
    old_rpkms = rpkms * (1-pig_matrix)

    new_counts, old_counts, summed_tot_counts = fill_true_zeros(counts, new_counts, old_counts, rounded = True)
    new_rpkms, old_rpkms, summed_tot_rpkms = fill_true_zeros(rpkms, new_rpkms, old_rpkms, rounded = False)
    
    if outfolder:
        if not os.path.exists(os.path.join(outfolder)):
            os.mkdir(os.path.join(outfolder))
        
        print(f"Saving expression matrices to {outfolder}")
        new_counts.to_csv(os.path.join(outfolder, 'new_counts.txt'), sep = '\t')
        old_counts.to_csv(os.path.join(outfolder, 'old_counts.txt'), sep = '\t')
        summed_tot_counts.to_csv(os.path.join(outfolder, 'tot_counts.txt'), sep = '\t')

        new_rpkms.to_csv(os.path.join(outfolder, 'new_rpkms.txt'), sep = '\t')
        old_rpkms.to_csv(os.path.join(outfolder, 'old_rpkms.txt'), sep = '\t')
        summed_tot_rpkms.to_csv(os.path.join(outfolder, 'tot_rpkms.txt'), sep = '\t')
    
    return new_counts, old_counts, summed_tot_counts, new_rpkms, old_rpkms, summed_tot_rpkms

def fill_true_zeros(og_counts, new, old, rounded = False):
    common_genes = np.intersect1d(og_counts.index.values, new.index.values)
    common_conds = np.intersect1d(og_counts.columns.values, new.columns.values)
    subset_og = og_counts.loc[common_genes, common_conds]
    
    num_cols = subset_og.shape[1]

    zero_dict = {}
    for colidx in range(num_cols):
        rowindeces = np.where(subset_og.iloc[:, colidx] == 0)[0].tolist()        
        colname    = subset_og.columns[colidx]
        genenames  = subset_og.index[rowindeces].values.tolist()
        zero_dict[colname] = genenames
        
    for bc, gnames in zero_dict.items():
        new.loc[gnames, bc] = 0
        old.loc[gnames, bc] = 0
    
    tot = new + old
        
    if rounded:
        new = new.round(0)
        old = old.round(0)
        tot = tot.round(0)
        
    return new, old, tot

### Plotting functions ###

def plot_conversion_rates(conv_rates, outpath = None):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (15, 5))
    
    flierprops = dict(marker='o', markerfacecolor='None', markersize = 6,  markeredgecolor='black')
    
    sns.boxplot(x = 'conversionType', y = 'conversionRate', data = conv_rates.query("strand == '+'"),
                ax = ax1, flierprops = flierprops)
    sns.boxplot(x = 'conversionType', y = 'conversionRate', data = conv_rates.query("strand == '-'"),
                ax = ax2, flierprops = flierprops)
    
    [ax.tick_params(labelsize = 12) for ax in (ax1, ax2)]
    ax1.set_title('Positive strand genes', fontsize = 12)
    ax2.set_title('Negative strand genes', fontsize = 12)
    
    ax1.set_xlabel("Conversion type", fontsize = 12)
    ax2.set_xlabel("Conversion type", fontsize = 12)
    
    ax1.set_ylabel("Conversion rate", fontsize = 12)
    ax2.set_ylabel("Conversion rate", fontsize = 12)
    
    fig.suptitle("Conversion rates", fontsize = 15)
    sns.despine()
    
    if outpath:
        if not os.path.exists(os.path.split(outpath)[0]):
            os.makedirs(os.path.split(outpath)[0])
        plt.savefig(outpath)

def plot_pig_read_depth_dependence(pigs, frac = 0.001, cmap_name = "RdBu_r", outpath = None):
    sample = pigs.sample(frac = frac)

    fig, ax = plt.subplots(figsize = (15, 7))

    x  = sample.numreads
    y  = sample.conf_width
    xy = np.vstack([x,y])
    z  = gaussian_kde(xy)(xy)
    
    ax.scatter(x, y, c = z, alpha = 0.1, cmap = cmap_name)
    ax.set_ylim([0, 1])
    
    norm = mpl.colors.Normalize(vmin=min(z), vmax=max(z))
    
    divider = make_axes_locatable(ax)
    cax = divider.new_vertical(size='5%', pad=0.6, pack_start = True)
    fig.add_axes(cax)
    
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap_name),
                 cax=cax, orientation='horizontal', label='Normalized density')
    
    ax.set_xlabel('Number of reads', fontsize = 12)
    ax.set_ylabel('$\pi_g$ estimate confidence interval', fontsize = 12)
    
    ax.tick_params(labelsize = 12)
    ax.set_xscale('log')
    ax.set_title(f"Relationship between number of reads sequenced and $\pi_g$ confidence interval width\nPlotted using {frac*100}% of the data")
    sns.despine()
    
    if outpath:
        print(f"Saving file to {outpath}...")
        if not os.path.exists(os.path.split(outpath)[0]):
            os.makedirs(os.path.split(outpath)[0])
        plt.savefig(outpath, bbox_inches = "tight")
            
    plt.show()

def plot_readstats(stats_df, readcount_threshold = 200_000, ngene_threshold = 7000, exon_threshold = 45, snr_threshold = 10, ngene_accessor = 'ngenes_with_pigs', outpath = None):
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (25, 5))

    ax1.scatter(stats_df.total_reads, (stats_df.Exon / stats_df.total_reads) * 100,
              s = 100, ec = 'k', color = 'C0', alpha = 0.5)
    ax1.tick_params(labelsize = 12)
    
    ax1.axvline(readcount_threshold, linestyle = '--', c = 'k', alpha = 0.3)
    ax1.axhline(exon_threshold, linestyle = '--', c = 'k', alpha = 0.3)
    
    ax1.set_xlabel('Total number of reads', fontsize = 14)
    ax1.set_ylabel('Percent exonic reads', fontsize = 14)
    ax1.set_xscale("log")
    
    if not ngene_accessor in ["ngenes_with_pigs", "ngenes_with_pigs_gtzero"]:
        raise ValueError("Error: ngene_accessor needs to be one of [ngenes_with_pigs/ ngenes_with_pigs_gtzero]")
    
    ax2.scatter(stats_df.total_reads,
                stats_df.loc[:, ngene_accessor],
                s = 110, ec = 'k', alpha = 0.7, c = 'coral')
    
    ax2.set_xlabel('Total number of reads', fontsize = 14)
    ax2.set_ylabel('Number of genes with inferred $\pi_g$', fontsize = 14)
    
    ax2.tick_params(labelsize = 12)
    ax2.axvline(readcount_threshold, linestyle = '--', c = 'k', alpha = 0.5)
    ax2.axhline(ngene_threshold, linestyle = '--', c = 'k', alpha = 0.5)
    ax2.set_xscale('log')
    
    ax3.scatter(stats_df.total_reads, stats_df.snr,  ec = 'k', s = 110, alpha = 0.7, c = 'mediumorchid')
    ax3.tick_params(labelsize = 12)
    ax3.set_xlabel('Total number of reads', fontsize = 12)
    ax3.set_ylabel('Signal to noise ratio', fontsize = 12)
    
    ax3.axhline(snr_threshold, linestyle = '--', c = 'k', alpha = 0.5)
    ax3.axvline(readcount_threshold, linestyle = '--', c = 'k', alpha = 0.5)
    ax3.set_xscale('log')
    
    fig.suptitle('Condition filtering', fontsize = 18)

    out_dict = {"readcount_threshold": readcount_threshold,
                "ngene_threshold": ngene_threshold,
                "exon_threshold": exon_threshold,
                "snr_threshold": snr_threshold,
                "ngene_accessor": ngene_accessor}
    
    sns.despine()
    
    if outpath:
        print(f"Saving file to {outpath}...")
        if not os.path.exists(os.path.split(outpath)[0]):
            os.makedirs(os.path.split(outpath)[0])
        plt.savefig(outpath)
        
    plt.show()
    return out_dict
