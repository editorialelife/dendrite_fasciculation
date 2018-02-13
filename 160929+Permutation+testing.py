
# coding: utf-8

# # Permutation testing (for genotypes compared to WT)

# In[1]:
import numpy
import os
from csv import reader
from numpy import random
import pylab

def permutation_test(num_of_animals, WT, mutant, timepoint, number_of_loops):
    """
    input: two lists - WT and genotype - distance * # animals
    returns: plot of percentiles v. distance
    
    1. Merge populations between genotype and WT data (needs to be equal number of bins)
    
    for i in range(500-1000):
    2. Make WT and genotype populations (index=0)
    3. Random sampling without replacement to make two populations
    4. Get count_ones, count_twos, count_threes, total_count for (a) experimental and (b) expected
    5. run fisher exact test --> get p-values for simulated and true population sets
    Output will be a matrix where rows = pixels and columns = number of comparisons
    
    6. for each pixel: use scipy.stats.percentile of score - percentileofscore(list, score) and that returns 
    the percentile the 'true' comparison is in
    7. graph the percentiles - x-axis is distance and y-axis is percentile

    """
    
    # get control and mutant genotypes middle dendrite points for population
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    
    # generate a random dataset
    def random_dataset(n=21):
        random_sets = []
        for i in range(n):
            random_set = []
            for i in range(98):
                random_set.append(str(random.randint(1,4)))
            random_sets.append(random_set)            
        return random_sets
    
        
        
        
    if WT == 'random':
        control_animals = random_dataset(num_of_animals)
    else:
        rel_path_WT = 'SADA csv files binned/' + WT + '/' + WT + '_' + timepoint + '_middle_dendr_pop.csv'
        abs_file_path_WT = os.path.join(script_dir, rel_path_WT)
        control = open(abs_file_path_WT, 'r')
        
        csv_control = reader(control)
        next(csv_control)
        control_animals = []
        for control_animal in control:
            control_animal = control_animal[:-3]
            control_animal = control_animal.split(',')
            control_animal = control_animal[1:]
            
#            # for WT data, your 123 key is different than other data, remember?
#            if WT == 'WT':
#                mydict = {'0':'0', '1':'3', '2':'1', '3':'2'}
#                for i in range(len(control_animal)):
#                    control_animal[i] = mydict[control_animal[i]]
            control_animals.append(control_animal)
        control.close()
    control_animals = control_animals[:num_of_animals]

        
    
    if mutant == 'random':
        experimental_animals = random_dataset(len(control_animals))
    else:
        rel_path_exp = 'SADA csv files binned/' + mutant + '/' + mutant + '_' + timepoint + '_middle_dendr_pop.csv'
        abs_file_path_exp = os.path.join(script_dir, rel_path_exp)
        experimental = open(abs_file_path_exp, 'r')
        
        csv_experimental = reader(experimental)
        next(csv_experimental)
        experimental_animals = []
        for experimental_animal in experimental:
            experimental_animal = experimental_animal[:-3]
            experimental_animal = experimental_animal.split(',')
            experimental_animal = experimental_animal[1:]
            experimental_animals.append(experimental_animal)
        experimental.close()
    experimental_animals = experimental_animals[:num_of_animals]
    
    
    # merge control and mutant genotypes
    mixed_genotypes = control_animals + experimental_animals
    assert (len(control_animals) + len(experimental_animals) == len(mixed_genotypes)), 'something went wrong with adding the genotypes'
    
    print(len(control_animals))
    print(len(experimental_animals))
    print(len(mixed_genotypes))
    
    # get fisher test p values for true control and experimental populations
    fisher_test_true_population = fisher_test(control_animals, experimental_animals)
    
    # loop through number_of_loops times for permutation testing
    fisher_test_populations = []
    for i in range(number_of_loops):
        # Make fake WT population, fake_control
        fake_experimental = mixed_genotypes[:]
        fake_control_index = random.choice(len(mixed_genotypes), len(control_animals), replace=False)
        fake_control = []
        for number in fake_control_index:
            fake_control.append(fake_experimental[number])
        
        # Make fake experimental population, fake_experimental
        fake_control_index_sorted = sorted(fake_control_index, reverse=True)
        for fake_control_index_sorted_number in fake_control_index_sorted:
            del(fake_experimental[fake_control_index_sorted_number])
        
        #get fisher test p values
        fisher_test_populations.append(fisher_test(fake_control, fake_experimental))
            
    # fisher_test_population is population data including "fake" fisher tests
    fisher_test_populations.insert(0,fisher_test_true_population) # [0] is the true population p values

    import pandas as pd
    fisher_df = pd.DataFrame(fisher_test_true_population)
    fisher_df.to_csv('WT_fisher_pvals.csv')

    # get distributions graphs
    get_normal_distributions(WT, mutant, timepoint, number_of_loops, fisher_test_populations)

    # get percentile score for true population
    fisher_percentile = get_percentile(fisher_test_populations)
    
    # get smoothed values
    fisher_percentile_smoothed = smoothing(fisher_percentile, 5)
    
    # save percentile arrays to csv
    import pandas as pd
    
    df_fisher_smooth = pd.DataFrame(fisher_percentile_smoothed)
    df_fisher = pd.DataFrame(fisher_percentile)
    df_fisher_smooth.to_csv('./SADA csv files binned/' + WT + '_' + mutant + '_' + timepoint + '_' + str(number_of_loops) + '_fisherpercentiles_smoothed.csv')
    df_fisher.to_csv('./SADA csv files binned/' + WT + '_' + mutant + '_' + timepoint + '_' + str(number_of_loops) + '_fisherpercentiles_values.csv')

    # get percentile_plot and bar
    fisher_percentile_smoothed
    percentile_plot(fisher_percentile_smoothed, WT, mutant, timepoint, number_of_loops)
    percentile_bar(fisher_percentile_smoothed, mutant)
    pylab.xlabel('distance along dendrite')
    pylab.ylabel('percentile rank')
    pylab.title('Fisher test p-value ranks')
    pylab.savefig('./SADA csv files binned/' + WT + '_' + mutant + '_' + timepoint + '_' + str(number_of_loops) + '_fisherpercentiles.pdf')
    pylab.close()
    
    return fisher_percentile
    
    
def fisher_test(control, experimental):
    """
    input: two lists, control and experimental middle dendrite values for two populations, along the length of the dendrite
    returns: fisher test p-values
    """
    # get scaling factor, control/experimental, will multiply to experimental count
    scaling_factor = len(control)/len(experimental)
    
    # control and experimental have 19 or 27 rows and 100 columns; when transposed they will have 100 rows x 19 or 27 cols
    control_transposed = [list(x) for x in zip(*control)]
    experimental_transposed = [list(x) for x in zip(*experimental)]

    # get counts for '1', '2', '3', '0' for experimental and controls
    control_ones = []
    control_twos = []
    control_threes = []
    control_zeroes = []
    for control_bin in control_transposed:
        control_one = control_bin.count('1')
        control_ones.append(control_one)
        control_two = control_bin.count('2')
        control_twos.append(control_two)
        control_three = control_bin.count('3')
        control_threes.append(control_three)
        control_zero = control_bin.count('0')
        control_zeroes.append(control_zero)
        
    control_values = [control_ones, control_twos, control_threes, control_zeroes] # count of ones, twos, threes, zeroes at each pixel
    
    
    experimental_ones = []
    experimental_twos = []
    experimental_threes = []
    experimental_zeroes = []
    for experimental_bin in experimental_transposed:
        experimental_one = experimental_bin.count('1') * scaling_factor
        experimental_ones.append(experimental_one)
        experimental_two = experimental_bin.count('2') * scaling_factor
        experimental_twos.append(experimental_two)
        experimental_three = experimental_bin.count('3') * scaling_factor
        experimental_threes.append(experimental_three)
        experimental_zero = experimental_bin.count('0') * scaling_factor
        experimental_zeroes.append(experimental_zero)

    experimental_values = [experimental_ones, experimental_twos, experimental_threes, experimental_zeroes] 

    
    # check your scaling factor was ok
    total_control_check = control_ones[0] + control_twos[0] + control_threes[0] + control_zeroes[0]
    total_experimental_check = experimental_ones[0] + experimental_twos[0] + experimental_threes[0] + experimental_zeroes[0]
    assert (abs(total_experimental_check - total_control_check) < 1), 'something is wrong with your scaling factor'

    
    # fisher 2x3:
    import matlab.engine
    eng = matlab.engine.start_matlab()
    try:
        eng.addpath(r'/Users/ZhiqiYip/Documents/MATLAB/Code')
    except AttributeError:
        eng.addpath(r'/Users/hink_pink/Documents/MATLAB/Code')

    # make MATLAB array
    fisher_p_values = []
    for i in range(len(experimental_values[0])):
        
        array_for_fisher = matlab.double([
                [int(control_values[0][i]), int(control_values[1][i]), int(control_values[2][i])],
                [int(experimental_values[0][i]), int(experimental_values[1][i]), int(experimental_values[2][i])]
            ])
        
        # get fisher p-values
        fisher_p_value = eng.myfisher23(array_for_fisher)
        fisher_p_values.append(fisher_p_value)   

    return fisher_p_values



def get_percentile(fisher_values_compiled):
    """
    input: fisher_values_compiled, where [0] is the p values of the true population and [1:] is the p values of the "fake" populations
    returns: percentiles of true population values
    """
    from scipy import stats
    percentiles = []
    for i in range(len(fisher_values_compiled[0])):
        p_values_by_column = []
        
        for population in fisher_values_compiled:
            p_values_by_column.append(population[i])

        percentile = stats.percentileofscore(p_values_by_column, p_values_by_column[0])
        percentiles.append(percentile)
    
    return percentiles


    
    
def get_normal_distributions(control, experimental, timepoint, number_of_loops, fisher_values_compiled):
    """
    input: fisher_values_compiled, where [0] is the p values of the true population and [1:] is the p values of the "fake" populations
    returns: plots showing distribution of p-values for each point along the dendrite and where the "true" p-value is located.
    """
    dendrite_locations_to_sample = list(range(5, 100, 10))
    
    # get distribution at many locations
    for location in dendrite_locations_to_sample:
        distr_values = []
        for sample in fisher_values_compiled:
            distr_values.append(sample[location])
        distr_values = numpy.array(distr_values)
        
        # make actual figure
        pylab.figure()
        pylab.hist(distr_values, 10)
        pylab.plot(distr_values[0], 1, 'ro')
        
        # reverse the axes
        ax = pylab.gca()
        ax.set_xlim(ax.get_xlim()[::-1])
                
        pylab.title('distribution of values at bin ' + str(location))
        pylab.xlabel('fisher p-values')
        pylab.ylabel('count')
        pylab.savefig('./SADA csv files binned/' + control + '_' + experimental + '_' + timepoint + '_' + str(number_of_loops) + '_distribution_at_' + str(location) + '.pdf')
        pylab.close()
    
    
def percentile_plot(percentiles, control, experimental, timepoint, number_of_loops):
    """
    input: fisher percentiles of true values
    returns: None
    prints: percentile plot of true values
    """
    
    # make percentile plot
    pylab.figure()
    pylab.plot(percentiles)


def percentile_bar(percentiles, experimental):
    """
    input: fisher percentiles of true values
    prints: percentile histogram
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
        
    # make percentile colorbar
    intensity = numpy.array([percentiles]) # this nested list is important
    
    # cmap
    if experimental == 'random':
        cmap = 'binary_r' # fisher test against random
    else:
        cmap = 'binary_r' # fisher test against genotype

    plt.pcolor(intensity, norm=LogNorm(vmin=1/501*100, vmax=100.0), cmap=cmap)
    plt.colorbar()


def smoothing(values, window):
    """
    Compute the moving average of y with specified window length.
    Args:
        y: an 1-d pylab array with length N, representing the y-coordinates of the N sample points
        window_length: an integer indicating the window length for computing moving average
    Returns:
        a list with the same length as y storing moving average of y-coordinates of the N sample points
    """
    # create list of smoothed data points
    averages = []
    
    for i in range(len(values)):
        
        # create list of data contained in window
        last_y_slice = i + window + 1
        if last_y_slice >= len(values) + 1:
            last_y_slice = len(values) + 1
        
        sma = pylab.array(values[i:last_y_slice])
        
        # get average of data in window
        average = sma.mean()
        averages.append(average)
    
    # return array of moving averages
    return averages    
    
    
#%%
# permutation_test(9, 'ok1489', 'WT', '48h', 500)
permutation_test(20, 'WT', 'random', '48h', 500)

# permutation_test(15, 'ok711', 'WT', '48h', 500)
# permutation_test(15, 'ok711', 'random', '48h', 500)
# tests
## random tests
# random_list = ['e1745', 'eq1', 'hmnIs23 M153L2', 'ky146', 'ky146mu256', 'm86sa204', 'M132L4', 'M144L2 N2', 
#                 'M157L1', 'mu256', 'ok244', 'ok711', 'rh310', 'rh310_control', 'u74', 'WT', 'wy686']

# for random_l in random_list:
#     permutation_test(random_l, 'random', '48h', 500)

# timepoints = ['WT', 'ky146', 'mu256', 'rh310']
# hours = ['24h', '72h']
# for timepoint in timepoints:
#     for hour in hours:
#         permutation_test(timepoint, 'random', hour, 500)

# print('random tests done')


# ## WT tests
# WT_list = ['e1745', 'eq1', 'hmnIs23 M153L2', 'ky146', 'ky146mu256', 'm86sa204', 'M132L4', 'M144L2 N2', 
#                 'M157L1', 'mu256', 'ok244', 'ok711', 'rh310', 'rh310_control', 'u74', 'wy686']
    
# for WT_l in WT_list:
#     permutation_test(WT_l, 'WT', '48h', 500)

# WT_timepoints = ['ky146', 'mu256', 'rh310']
# for WT_timepoint in WT_timepoints:
#     for hour in hours:
#         permutation_test(timepoint, 'WT', hour, 500)

# print('WT tests done')


# ## other weird ones
# ### double mutant
# permutation_test('ky146mu256', 'ky146', '48h', 500)
# permutation_test('ky146mu256', 'mu256', '48h', 500)

# ### rh310_control
# permutation_test('rh310_control', 'rh310', '48h', 500)

# print('other tests done')


print ('all done!')