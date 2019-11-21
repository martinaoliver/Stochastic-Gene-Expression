import matplotlib.pyplot as plt
import statistics
from gillespiefunctions import celldivision_gillespie_total
import random
import numpy as np

##Run Gillespie with Binomial Partitioning
#dt,dm,dd are dictionaries for storing gillespie outputs (time and mRNA). Every simulation will be stored in an item of this dictionary
dt={}
dm={}
dd = {}
simulation_number = 200
time_range = 10000
for n in range (1,simulation_number+1):
 dt['t_%s' % n], dm['m_%s' % n], dd['md_%s' % n] = celldivision_gillespie_total(time_range)
#%%
#plot
plt.figure(1)
plt.plot(dt['t_1'],dd['md_1'], 'g', linewidth=0.9)
plt.xlabel('Time (s)', size=10)
plt.ylabel('mRNA units', size=10)
plt.savefig('4,celldiv_tot')
plt.show()
plt.close()

#%%
##Statistical Analysis
#fano function performs a statistical analysis of the Gillespie data. The mean and variance is calculated accross the simulation
#Then, all means and variances are amongst all simulations
#The fano factor is calculated from the averaged mean and variance. 
def fano(dictionary, dataset):
    mean_l = []
    variance_l = []
    for n in range (1,simulation_number+1):
        mean = statistics.mean(dictionary[dataset + '_%s' %n][100:-1])
        variance = statistics.variance(dictionary[dataset + '_%s' %n][100:-1])
        mean_l.append(mean)
        variance_l.append(variance)
    mean_simulations = statistics.mean(mean_l)
    variance_simulations = statistics.mean(variance_l)
    fano_simulations = variance_simulations/mean_simulations
    return (mean_simulations, variance_simulations, fano_simulations, mean_l)
    print (mean_simulations, variance_simulations, fano_simulations)

mean_simulations, variance_simulations, fano_simulations, mean_l = fano(dd,'md')
print (mean_simulations, variance_simulations, fano_simulations)


#%%
plt.figure(2)
plt.hist(mean_l, 40)
plt.xlabel('mRNA values')
plt.ylabel('Frequency')
plt.savefig('4,mean_tot.png')
plt.close()
