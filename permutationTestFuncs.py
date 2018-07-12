import numpy as np
import os
from csv import reader
from numpy import random
from matplotlib import pyplot as plt
import pandas as pd

# Functions for permutationTest

def saveFisherCsvAndPlots(fisher_percentile, fisher_percentile_smoothed, WT, mutant, timepoint, number_of_loops):
    df_fisher = pd.DataFrame(fisher_percentile, columns=['fisher_percentile'])
    df_fisher['fisher_smoothed_percentile'] = fisher_percentile_smoothed
    df_fisher.to_csv('./SADA csv files binned/' + WT + '_' + \
                     mutant + '_' + timepoint + '_' + \
                     str(number_of_loops) + \
                     '_fisherpercentiles_values.csv')

    # get percentile_plot and bar
    # fisher_percentile_smoothed
    percentile_bar(fisher_percentile_smoothed, mutant, WT)
    plt.xlabel('distance along dendrite')
    plt.ylabel('percentile rank')
    plt.title(WT + '_' + mutant + '_' + timepoint + '_' + \
                str(number_of_loops) + '_Fisher test p-value ranks')
    plt.savefig('./SADA csv files binned/' + WT + '_' + mutant + \
                  '_' + timepoint + '_' + str(number_of_loops) + \
                  '_fisherpercentiles.pdf')
    plt.close()
    
def random_dataset(n=50):
    """
    generate a random dataset
    """
    random_sets = []
    for i in range(n):
        random_set = []
        for i in range(98):
            random_set.append(random.randint(1,4))
        random_sets.append(random_set)  
    random_sets = pd.DataFrame(random_sets).T
    
    return random_sets

def get_percentile(fisher_values_compiled):
    """
    input: fisher_values_compiled, where [0] is the p values 
    of the true population and [1:] is the p values of the "fake" populations
    returns: percentiles of true population values
    """
    from scipy import stats
    percentiles = []
    df = pd.DataFrame(fisher_values_compiled, columns=None)
    for i in df.columns:
        percentile = stats.percentileofscore(df[i].values, df[i].values[0]) / 100.
        percentiles.append(percentile)

    # for i in range(len(fisher_values_compiled[0])):
    #     p_values_by_column = []
        
    #     for population in fisher_values_compiled:
    #         p_values_by_column.append(population[i])

    #     percentile = stats.percentileofscore(p_values_by_column, p_values_by_column[0])
    #     percentiles.append(percentile)
    
    return percentiles

def get_normal_distributions(control, experimental, timepoint, number_of_loops, fisher_values_compiled):
    """
    input: fisher_values_compiled, where [0] is the p values of the
    true population and [1:] is the p values of the "fake" populations
    returns: plots showing distribution of p-values for each point 
    along the dendrite and where the "true" p-value is located.
    """
    dendrite_locations_to_sample = list(range(5, 100, 10))
    
    # get distribution at many locations
    for location in dendrite_locations_to_sample:
        distr_values = []
        for sample in fisher_values_compiled:
            distr_values.append(sample[location])
        distr_values = np.array(distr_values)
        
        # make actual figure
        plt.figure()
        plt.hist(distr_values, 10)
        plt.plot(distr_values[0], 1, 'ro')
        
        # reverse the axes
        ax = plt.gca()
        ax.set_xlim(ax.get_xlim()[::-1])
                
        plt.title('distribution of values at bin ' + str(location))
        plt.xlabel('fisher p-values')
        plt.ylabel('count')
        plt.savefig('./SADA csv files binned/' + control \
                      + '_' + experimental + '_' + timepoint + \
                      '_' + str(number_of_loops) + '_distribution_at_'\
                      + str(location) + '.pdf')
        plt.close()

def smoothing(values, window):
    """
    Compute the moving average of y with specified window length.
    Args:
        y: an 1-d pylab array with length N, representing the 
        y-coordinates of the N sample points
        window_length: an integer indicating the window length for 
        computing moving average
    Returns:
        a list with the same length as y storing moving 
        average of y-coordinates of the N sample points
    """
    # create list of smoothed data points
    averages = []
    
    for i in range(len(values)):
        
        # create list of data contained in window
        last_y_slice = i + window + 1
        if last_y_slice >= len(values) + 1:
            last_y_slice = len(values) + 1
        
        sma = np.array(values[i:last_y_slice])
        
        # get average of data in window
        average = sma.mean()
        averages.append(average)
    
    # return array of moving averages
    return averages    
    
# Functions for saveFisherCsvAndPlots
def percentile_plot(percentiles, control, experimental, timepoint, number_of_loops):
    """
    input: fisher percentiles of true values
    returns: None
    prints: percentile plot of true values
    """
    # make percentile plot
    plt.figure()
    plt.plot(percentiles)
    plt.close()

def percentile_bar(percentiles, experimental, WT):
    """
    input: fisher percentiles of true values
    prints: percentile histogram
    """
    from matplotlib.colors import LogNorm
        
    # make percentile colorbar
    intensity = np.array([percentiles]) # this nested list is important
    
    # cmap
    if experimental == 'random':
        cmap = 'Blues_r' # fisher test against random
    elif WT=='random':
        cmap='Blues_r'
    else:
        cmap = 'Reds_r' # fisher test against genotype

    plt.pcolor(intensity, norm=LogNorm(vmin=1/501, vmax=1.), cmap=cmap)
    ax = plt.gca()
    ax.invert_yaxis()
    plt.colorbar()
    
def permutationTest(WT, mutant, timepoint, number_of_loops):
    if WT == 'random':
        control_animals = random_dataset()
    else:
        control_animals = pd.read_csv('SADA csv files binned/' + WT + '/' + \
                              WT + '_' + timepoint + '_middle_dendr_pop.csv',
                              index_col='Unnamed: 0').iloc[:98]

    if mutant == 'random':
        experimental_animals = random_dataset(control_animals.shape[1])
    else:
        experimental_animals = pd.read_csv('SADA csv files binned/' + mutant + '/' + \
                              mutant + '_' + timepoint + '_middle_dendr_pop.csv',
                              index_col='Unnamed: 0').iloc[:98]

    sample_size = min(control_animals.shape[1], experimental_animals.shape[1])

    control_animals = control_animals.replace(np.nan, 0).iloc[:, :sample_size]
    experimental_animals = experimental_animals.replace(np.nan, 0).iloc[:, :sample_size]

    # merge control and mutant genotypes
    mixed_genotypes = experimental_animals.join(control_animals) # each column is a different animal
    
    # get fisher test p values for true control and experimental populations
    fisher_test_true_population = fisher_test(control_animals, experimental_animals)

    # loop through number_of_loops times for permutation testing
    fisher_test_populations = []

    for i in range(number_of_loops):
        # Make fake control population
        fake_experimental = mixed_genotypes[:]
        fake_control_index = list(random.choice(mixed_genotypes.columns.values, 
                                           control_animals.shape[1], 
                                           replace=False))
        fake_control = mixed_genotypes[fake_control_index]

        # Make fake experimental population
        fake_experimental = mixed_genotypes.drop(fake_control_index, axis=1)

        # get fisher test p values
        fisher_test_populations.append(fisher_test(fake_control, fake_experimental))

    # fisher_test_population is population data including "fake" fisher tests
    fisher_test_populations.insert(0, fisher_test_true_population)

    # save fisher test p-values for true population
    fisher_df = pd.DataFrame(fisher_test_true_population, columns=['p_value'])
    fisher_df.to_csv(WT + '_' + mutant + '_' + timepoint + '_' + \
                     str(number_of_loops) + '_' + 'true_population_fisher_pvals.csv')

    # get distributions graphs
    get_normal_distributions(WT, mutant, timepoint, number_of_loops, fisher_test_populations)

    # get percentile score for true population
    fisher_percentile = get_percentile(fisher_test_populations)

    # get smoothed values
    fisher_percentile_smoothed = smoothing(fisher_percentile, 5)
    
    return fisher_percentile, fisher_percentile_smoothed, WT, mutant, timepoint, number_of_loops



def fisher_test(control, experimental):
    """
    input: two lists, control and experimental middle dendrite values 
    for two populations, along the length of the dendrite
    returns: fisher test p-values
    """
    control = control.astype(int)
    experimental = experimental.astype(int)
    
    # get counts for '1', '2', '3', '0' for experimental and controls at each pixel
    control_counts = control.apply(pd.value_counts, axis=1)
    experimental_counts = experimental.apply(pd.value_counts, axis=1)
    control_counts = control_counts.replace(np.nan, 0)
    experimental_counts = experimental_counts.replace(np.nan, 0)
    control_counts = control_counts.astype(int)
    experimental_counts = experimental_counts.astype(int)

    if control_counts.shape[1] > 3:
        control_counts = control_counts[[1, 2, 3]]
    if experimental_counts.shape[1] > 3:
        experimental_counts = experimental_counts[[1, 2, 3]]

    # fisher 2x3:
    import matlab.engine
    eng = matlab.engine.start_matlab()
    try:
        eng.addpath(r'/Users/ZhiqiYip/Documents/MATLAB/Code')
    except AttributeError:
        eng.addpath(r'/Users/hink_pink/Documents/MATLAB/Code')

    # make MATLAB array
    fisher_p_values = []
    
    for i in range(experimental_counts.shape[0]):
        control_row = control_counts.iloc[i].values
        experimental_row = experimental_counts.iloc[i].values
        
        array_for_fisher = matlab.double([
                            [int(control_row[0]), 
                             int(control_row[1]), 
                             int(control_row[2])],
                            [int(experimental_row[0]), 
                             int(experimental_row[1]), 
                             int(experimental_row[2])]
                            ])
        
        # get fisher p-values
        fisher_p_value = eng.myfisher23(array_for_fisher)
        fisher_p_values.append(fisher_p_value)   

    return fisher_p_values