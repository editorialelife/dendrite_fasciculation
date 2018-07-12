import pandas as pd
import numpy as np
from itertools import zip_longest
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# maximum pairwise distance csv + plot + middledendrite csv
def maxpairwiseAndMiddleDendrite(genotype, timepoint):
    max_pairwise_df = pd.DataFrame()
    middle_dendrite_df = pd.DataFrame()

    all_animals = pd.read_csv('180531_defasc_classifier_ver3.csv', 
                              index_col=['Unnamed: 0'])\
                    [['genotype','timepoint','id','test5']]

    for animal in all_animals[(all_animals.genotype==genotype) & (all_animals.timepoint==timepoint)].id.values:
        single_animal_df = pd.read_csv('./' + animal + '.csv')

        # cut off pixels toward cell body for all pairwise distances
        normalized_pairwise_distances = normalization(single_animal_df.RG_dist, 
                                                      single_animal_df.BR_dist, 
                                                      single_animal_df.GB_dist, 
                                                      animal, timepoint)
        binned_distances_transposed = pd.DataFrame(binned_dist(normalized_pairwise_distances[0], 
                                                   normalized_pairwise_distances[1], 
                                                   normalized_pairwise_distances[2], 
                                                   100), 
                                                   columns=['RG_binned', 'BR_binned', 'GB_binned'])
        binned_distances_transposed.to_csv('./' + animal + '_normalized_bins.csv')
        binned_distances_transposed['max_pairwise'] = binned_distances_transposed.max(axis=1)

        # get max pairwise distances
        max_pairwise_to_plot = all_animals[all_animals.id==animal].test5.tolist() + binned_distances_transposed['max_pairwise'].tolist()
        if len(all_animals[all_animals.id==animal].test5.tolist() + binned_distances_transposed['max_pairwise'].tolist()) == 100:
            max_pairwise_to_plot.append(0)
        max_pairwise_df[animal.split('/')[-1]] = max_pairwise_to_plot

        # get middle dendrite population.csv only if it's fasciculated
        if all_animals[all_animals.id==animal].test5.values=='bundled':
            switch_in_channel_order = {0:0, 1:3, 2:1, 3:2}
            binned_distances_transposed['middleDendrite'] = binned_distances_transposed.\
                                                            apply(middleDendrite, axis=1)

            # deal with switched channel order
            if animal.split('/')[-1].split('_')[1][:3] == '072':
                print(animal.split('/')[-1].split('_')[1] + ' changed channel order')
                binned_distances_transposed['middleDendrite'] = binned_distances_transposed.\
                                                                middleDendrite.\
                                                                apply(lambda x: switch_in_channel_order[x])
            elif animal.split('/')[-1].split('_')[1][:2] in ['05','06'] and animal.split('/')[-1].split('_')[1][4:6] == '18':
                print(animal.split('/')[-1].split('_')[1] + ' changed channel order')
                binned_distances_transposed['middleDendrite'] = binned_distances_transposed.\
                                                                middleDendrite.\
                                                                apply(lambda x: switch_in_channel_order[x])
            middle_dendrite_df[animal.split('/')[-1]] = binned_distances_transposed['middleDendrite']
        else:
            print(all_animals[all_animals.id==animal].id.values + 'is defasciculated')

    
    # max pairwise plots
    getMaxPairwisePlot(max_pairwise_df, timepoint=timepoint, genotype=genotype)
    defascThicknessPlot(animal=animal, max_pairwise_df=max_pairwise_df, genotype=genotype, timepoint=timepoint)
    
    max_pairwise_df.to_csv('./SADA csv files binned/' + genotype + '/' + \
                       genotype + '_' + timepoint + '_maxpairwise_pop.csv')
    middle_dendrite_df.to_csv('./SADA csv files binned/' + genotype + '/' + \
                       genotype + '_' + timepoint + '_middle_dendr_pop.csv')
    return middle_dendrite_df


def defascThicknessPlot(animal, max_pairwise_df, genotype, timepoint):
    animal = animal.lower()
    bundle_status = max_pairwise_df.iloc[0].values
    max_pairwise_df = max_pairwise_df.iloc[1:].reset_index(drop=True)

    if '60x' in animal:
        df_pos = max_pairwise_df/2/9.26
        df_neg = -max_pairwise_df/2/9.26
    elif '40x' in animal:
        df_pos = max_pairwise_df/2/6.646
        df_neg = -max_pairwise_df/2/6.646


    idx=0
    for i in max_pairwise_df.columns:
        if bundle_status[idx] == 'bundled':
            color='cornflowerblue'
        elif bundle_status[idx] == 'defasc':
            color='red'
        idx += 1
        plt.fill_between(df_pos.index.tolist(), df_pos.loc[:, i].tolist(), df_neg.loc[:, i].tolist(), 
                         color=color, 
                         alpha = 0.1,
                        )
    plt.ylim([-3, 3])
    plt.xlim([0, 100])
    plt.axhline(y=-1, color='black')
    plt.axhline(y=1, color='black')
    plt.grid(color='b', axis='y', alpha=0.5)
    plt.axis('off')
    plt.savefig('./SADA csv files binned/' + genotype + '/' + \
                       genotype + '_' + timepoint + '_defascThicknessPlot.svg')
    plt.close()

# heatmap
def getHeatmap(middle_dendrite_df, genotype, timepoint):
    from matplotlib.colors import LinearSegmentedColormap
    vmax = 3.0    
    cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / vmax, 'white'),
                                                    (1 / vmax, 'yellow'),
                                                    (2 / vmax, 'red'),
                                                    (3 / vmax, 'blue')])
    fig, ax = plt.subplots()
    plt.pcolor(middle_dendrite_df.T, cmap=cmap, vmin=0, vmax=vmax)
    ax.set_yticks(np.arange(middle_dendrite_df.shape[1]) + 0.5, minor=False)
    ax.invert_yaxis()
    ax.set_yticklabels(middle_dendrite_df.columns, minor = False)
    ax.tick_params('off')
    plt.title(genotype + '_' + timepoint + ' n=' + str(middle_dendrite_df.shape[1]))
    plt.savefig('./SADA csv files binned/' + genotype + '/' + \
                       genotype + '_' + timepoint + '_heatmap.svg')
    plt.close()

# make summary plot
def summaryPlot(middle_dendrite_df, genotype, timepoint):
    summary_plot_df = pd.DataFrame()
    for i in middle_dendrite_df.T.columns[::-1]:
        summary_plot_df = pd.DataFrame(middle_dendrite_df.T[i].value_counts()).join(summary_plot_df, how='outer')

    (summary_plot_df / summary_plot_df.sum()).T.\
    replace(np.nan, 0.).\
    rolling(5, center=True).mean().plot(color=['orange','r','b'])
    plt.legend('')
    plt.ylim([0.0,1.0])
    plt.xlim([0,100])
    plt.tick_params('off')
    plt.savefig('./SADA csv files binned/' + genotype + '/' + \
                       genotype + '_' + timepoint + '_summaryplot.svg')
    plt.close()
    
    return summary_plot_df

# chi-squared test + color bar
def chiSquareTest(summary_plot_df, genotype, timepoint):
    chi_sq_df = summary_plot_df.T
    chi_sq_df['total'] = chi_sq_df.sum(axis=1)
    chi_sq_df['expected_values'] = chi_sq_df.total / 3.
    try:
        chi_sq_df.columns=['zeroes', 'one','two','three','total','expected_values']
    except ValueError:
        chi_sq_df.columns=['one','two','three','total','expected_values']
    chi_sq_df = chi_sq_df.replace(np.nan, 0.)

    from scipy.stats import chisquare
    def getChiSquare(df):
        return chisquare([df['one'], df['two'], df['three']], df.expected_values)[1]

    chi_sq_df['chi_sq_p_vals'] = chi_sq_df.apply(getChiSquare, axis=1)
    chi_sq_df.rolling(5, center=True).mean().to_csv('./SADA csv files binned/' + genotype + '/' + \
                       genotype + '_' + timepoint + '_chisq.csv')

    plt.pcolor(pd.DataFrame(chi_sq_df.chi_sq_p_vals.\
                            rolling(5, center=True).mean().dropna()), 
               norm=LogNorm(vmin=0.000001, vmax=1.0), 
               cmap='Greens_r')
    plt.colorbar()
    plt.savefig('./SADA csv files binned/' + genotype + '/' + \
                       genotype + '_' + timepoint + '_chi_sq_bar.svg')
    plt.close()

# other fns
def middleDendrite(df):
    if df.RG_binned == df.BR_binned == df.GB_binned:
        return 0
    elif df.RG_binned == df.max_pairwise:
        return 1 # blue
    elif df.BR_binned == df.max_pairwise:
        return 2 # yellow
    elif df.GB_binned == df.max_pairwise:
        return 3 # red
    return 0

def middledendrite(binned_distances_transposed, file_name):
    """
    looking for middle dendrite points as defined by the point across the longest pairwise distance
    input: binned_distances_transposed
    returns: mid_dendrite_col
    """

    for row in binned_distances_transposed:
        if row[0] == row[1] == row[2]:
            mid_dendrite = 0
        elif row[0] == max(row):
            mid_dendrite = 1
        elif row[1] == max(row):
            mid_dendrite = 2
        elif row[2] == max(row):
            mid_dendrite = 3
        else: 
            mid_dendrite = 0

        # convert WT values that start with 072 to the right color order
        if file_name[3:6] == '072':
            mid_dendrite = mydict[mid_dendrite]                    

        # convert WT values that start with 072 to the right color order
        if file_name[3:5] == '05' and file_name[7:9] == '18':
            mid_dendrite = mydict[mid_dendrite]                    
        mid_dendrite_col.append(mid_dendrite)


    return mid_dendrite_col

def getMaxPairwisePlot(max_pairwise_df, timepoint, genotype):
    plt.figure()
    max_pairwise_df = max_pairwise_df.iloc[1:, :]
    max_pairwise_df.rolling(5, center=True).mean().plot(color='cornflowerblue')
    plt.legend('')

    if timepoint == '48h' or timepoint == '24h':
        plt.ylim([0.0, float(9*5)]) # 5 micron y-range at 60x objective
        plt.yticks(np.arange(0, float(9*5)+1, 9))

    elif timepoint == '72h' and genotype == 'WT':
        plt.ylim([0.0, float(6*5)]) # 5 micron y-range at 40x objective
        plt.yticks(np.arange(0, float(6*5)+1, 6))

    elif timepoint == '72h':
        plt.ylim([0.0, float(9*5)]) # 5 micron y-range at 60x objective
        plt.yticks(np.arange(0, float(9*5)+1, 9))

    else:
        print ('warning, y axis numbers for maxpairwise dist plot incorrect')
    plt.xlabel('distance along dendrite, 100 bins')
    plt.ylabel('maximum pairwise distance, 5 microns')
    plt.title('maxpairwise distance plot for ' + genotype + '_' + timepoint)

    plt.savefig('./SADA csv files binned/' + genotype + '/' +                 genotype + '_' + timepoint + '_maxpairwise_pop_plot.svg')
    plt.close()

def binned_dist(RG_dist, BR_dist, GB_dist, bins=100):
    pixels = len(RG_dist)
    bin_size = pixels/bins
    n=0
    bin_index = []
    while n <= pixels:
        bin_index.append(round(n))
        n += bin_size

    RG_binned = []
    BR_binned = []
    GB_binned = []

    n_bins = 0
    while n_bins < (len(bin_index)-1):
        RG_bin = sum(RG_dist[bin_index[n_bins] : bin_index[n_bins + 1]]) / (-bin_index[n_bins] + bin_index[n_bins + 1])
        RG_binned.append(RG_bin)
        BR_bin = sum(BR_dist[bin_index[n_bins] : bin_index[n_bins + 1]]) / (-bin_index[n_bins] + bin_index[n_bins + 1])
        BR_binned.append(BR_bin)
        GB_bin = sum(GB_dist[bin_index[n_bins] : bin_index[n_bins + 1]]) / (-bin_index[n_bins] + bin_index[n_bins + 1])
        GB_binned.append(GB_bin)
        n_bins += 1

    binned_distances = RG_binned, BR_binned, GB_binned
    binned_distances_transposed = np.transpose(binned_distances)

    return binned_distances_transposed

def normalization(RG_dist, BR_dist, GB_dist, animal, timepoint):
    normalization_factor_microns = int(animal.split('/')[-1].split('_')[0])

    # adjust normalization factor according to image magnification (40x/60x/100x)
    animal_lower = animal.lower()
    if '100x' in animal_lower:
        normalization_factor_pixels = 15 * normalization_factor_microns
    elif '40x' in animal_lower:
        normalization_factor_pixels = 6 * normalization_factor_microns
    elif '60x' in animal_lower and timepoint == '48h': # 60x
        normalization_factor_pixels = int(9 * normalization_factor_microns * 0.9) # 90% of the anatomical distance
    elif '60x' in animal_lower and timepoint == '24h':
        normalization_factor_pixels = int(9 * normalization_factor_microns * 0.8) # 80% of the anatomical distance                    
    elif '60x' in animal_lower and timepoint == '72h':
        normalization_factor_pixels = int(9 * normalization_factor_microns * 1.0) # 100% of the anatomical distance                    
    zeros = [0] * (normalization_factor_pixels)
    RG_dist_de = [x_rg + y_rg for x_rg, y_rg in zip_longest(zeros, RG_dist, fillvalue=0)]
    BR_dist_de = [x_br + y_br for x_br, y_br in zip_longest(zeros, BR_dist, fillvalue=0)]
    GB_dist_de = [x_gb + y_gb for x_gb, y_gb in zip_longest(zeros, GB_dist, fillvalue=0)]

    RG_dist_de = RG_dist_de[:normalization_factor_pixels]
    BR_dist_de = BR_dist_de[:normalization_factor_pixels]
    GB_dist_de = GB_dist_de[:normalization_factor_pixels]

    normalized_pairwise = RG_dist_de, BR_dist_de, GB_dist_de
    return normalized_pairwise

# Final
def applyAll(genotype, timepoint):
    middle_dendrite_df = maxpairwiseAndMiddleDendrite(genotype, timepoint)
    getHeatmap(middle_dendrite_df, genotype, timepoint)
    summary_plot_df = summaryPlot(middle_dendrite_df, genotype, timepoint)
    chiSquareTest(summary_plot_df, genotype, timepoint)