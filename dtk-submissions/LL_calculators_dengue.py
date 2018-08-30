import math

import numpy as np
import math
import random
from scipy.stats import binom, multinomial
from scipy.stats import dirichlet
from scipy.special import betaln, comb, gammaln


##GOAL: Calculate the log binomial coefficient for two integer inputs
##AUTHOR: K. James Soda
##DATE STARTED: June 30, 2017
##INPUT
#n = int; the total number of options to form a combination
#x = int; the number of elements in a combination
##OUTPUT
#float; the log binomial coefficient for the two integers 
##TESTED USING THE FOLLOWING CALLS:
#from LL_calculators_dengue import logCombination
#Raise TypeError exception
#logCombination(1.5,1)
#logCombination(5,2.9)
#Raise ValueError exception
#logCombination(-5,2)
#logCombination(5,-2)
#logCombination(2,5)
#Successful call
#logCombination(5,2)
##OPERATES TO EXPECTATION

def logCombination(n,x):
	#Check to make sure the inputs make sense
	if type(n) != int or type(x) != int:
		raise TypeError("logCombination can only handle int arguments")
	if n < 0 or x < 0:
		raise ValueError("logCombination can only handle non-negative arguments")
	if x > n:
		raise ValueError("In logCombination n must be larger than x")

	result = 0
	#Calculate the log of the numerator
	for i in range(0,x):
		result += math.log(n - i)

	#Subtract away the log of the denominator
	for j in range(1,x+1):
		result -= math.log(j)

	return result


##GOAL:  For a simulation arising from a known model, calculate the log-likelihood of the model based on a provided time 
##	series given that the time series originated under a binomial distribution whose probability of success stems from
##	a beta distribution.  In other words, the log-likelihood of the model based on a beta-binomial distribution.  The
##	original intention was to calculate a model's likelihood based on a time series of reported cases (i.e., empirical 
##	successes) provided the simulated population size (i.e., number of trials) and simulated new reported cases (i.e.,
##	simulated successes).
##AUTHOR: K. James Soda
##DATE STARTED: June 28, 2017
##INPUT
#raw_data = list of lists (int);  outer list should be a vector time series with time steps in some temporal unit (e.g., a year).  
#	Each element in the inner list (i.e., vector observation in the outer list's vector time series) should be a univariate
#	time series with step sizes smaller than the outer series (e.g., a month).  Each of these observations should be the number
#	of successes in the training data. The routine will concatenate the inner list into a single list in the order that the 
#	outer list specifies, so the key point is that this final list must form a coherent series.  For example, [[1,2,3],[4,5,6]] 
#	will be concatenated into [1,2,3,4,5,6].  Alternatively, provide a single univariate series as a list within a 
#	single-element list (e.g., [[1,2,3,4,5,6]]).
#sim_data = list of lists (int);  the two lists will have separate interpretations, but both are time series with over the same
#	period and with the same step size as the final, concatenated list from raw_data (see above).
#		sim_data[0] = list (int); each element is the number of trials for that time step.
#		sim_data[1] = list (int); each element is the number of simulated successes for that time step.
##OUTPUT
#float; the log-likelihood of the model based on the data in raw_data
##TESTED USING THE FOLLOWING CALLS:
#from LL_calculators_dengue import betaBinomial_dengue
#
#raw_data = [[3,4],[0,9],[2,3]]
#sim_data = [[10,12,11,15,12,11],[5,7,2,8,1,0]]
#betaBinomial_dengue(raw_data,sim_data)
##OPERATES TO EXPECTATION (see pg. 27-28 of the project notes, vol. 1 [purple]; the final likelihood should be -13.11532)

def betaBinomial_dengue_modified(raw_data, sim_data, pop_data):
	##Initial data processing
	#raw_data is provided as a list of lists if multiple years are included, but sim_data 
	#	is just a list.  Here we convert raw_data to a single list, if necessary.

    empiricalSeries = raw_data[0]
    simSeries = sim_data[0]
    popSeries = pop_data[0]

    if len(raw_data) > 1:
        for i in range(1,len(raw_data)):
            empiricalSeries.extend(raw_data[i])
            simSeries.extend(sim_data[i])
            popSeries.extend(pop_data[i])
    #print empiricalSeries

	##Check to make sure that raw_data and sim_data have the same length
    if len(empiricalSeries) != len(simSeries):
        raise RuntimeError("raw_data and sim_data[1] must have the same length to use betaBinomial_dengue.")

    ##Implement the log-likelihood calculations.  The assumption is that popSeries provides the
	##	the number of trials, simSeries provides the number of successes in the simulation, and raw_data
	##	provides the number of successes in the empirical data.  
    logLik = 0
    for i in range(0,len(empiricalSeries)):
        #Get the current alpha and beta hyperparameters
        n = int(popSeries[i])
        x_hat = int(simSeries[i])
        x = int(empiricalSeries[i])
        tempAlpha = 1 + x_hat
        #print "Alpha at " + str(i) + ": " + str(tempAlpha)
        tempBeta = 1 + n - x_hat
		#print "Beta at " + str(i) + ": " + str(tempBeta)

		#NOTE: betaln actually gives the naturual log of the beta function's absolute value, but this should
		#	not matter too much (I think) because the calls to beta should return positive values anyway.
		#	The comparable scenario arises in dirichlet_dengue with gammaln
		#print "logLik at " + str(i) + ": " + str(math.log(comb(n,x)) + betaln(x + tempAlpha,n - x + tempBeta) - betaln(tempAlpha,tempBeta)) + "\n\n"
        logLik += logCombination(n,x) + betaln(x + tempAlpha,n - x + tempBeta) - betaln(tempAlpha,tempBeta)

    return logLik

def betaBinomial_dengue_simple(raw_data, sim_data, pop_data):
	##Initial data processing
	#raw_data is provided as a list of lists if multiple years are included, but sim_data
	#	is just a list.  Here we convert raw_data to a single list, if necessary.

    empiricalSeries = raw_data
    simSeries = sim_data
    popSeries = pop_data

    if len(empiricalSeries) != len(simSeries):
        raise RuntimeError("raw_data and sim_data[1] must have the same length to use betaBinomial_dengue.")

    ##Implement the log-likelihood calculations.  The assumption is that popSeries provides the
	##	the number of trials, simSeries provides the number of successes in the simulation, and raw_data
	##	provides the number of successes in the empirical data.
    logLik = 0
    for i in range(0,len(empiricalSeries)):
        #Get the current alpha and beta hyperparameters
        n = int(popSeries[i])
        x_hat = int(simSeries[i])
        x = int(empiricalSeries[i])
        tempAlpha = 1 + x_hat
        #print "Alpha at " + str(i) + ": " + str(tempAlpha)
        tempBeta = 1 + n - x_hat
		#print "Beta at " + str(i) + ": " + str(tempBeta)

		#NOTE: betaln actually gives the naturual log of the beta function's absolute value, but this should
		#	not matter too much (I think) because the calls to beta should return positive values anyway.
		#	The comparable scenario arises in dirichlet_dengue with gammaln
		#print "logLik at " + str(i) + ": " + str(math.log(comb(n,x)) + betaln(x + tempAlpha,n - x + tempBeta) - betaln(tempAlpha,tempBeta)) + "\n\n"
        logLik += logCombination(n,x) + betaln(x + tempAlpha,n - x + tempBeta) - betaln(tempAlpha,tempBeta)

    return logLik



def multinomial_dengue(raw_data, sim_data):
    num_obs = len(raw_data)
    num_strains = len(raw_data[0])
    LL = 0.0

    for x in range(num_obs):
        temp_raw_data = list()
        temp_sim_data = list()
        for y in range(num_strains):
            if raw_data[x][y] > 0:
                temp_raw_data.append(raw_data[x][y])
                denominator = sum(sim_data[x][1:1+num_strains])
                if denominator == 0:
                    temp_sim_data.append(0)
                else:
                    temp_sim_data.append(sim_data[x][y+1]/denominator)
        if sum(temp_sim_data) == 0:
            length = len(temp_sim_data)
            proportion = 1. / length
            for y in range(length):
                temp_sim_data[y] = proportion
        LL_temp = multinomial.pmf(temp_raw_data, n=sum(temp_raw_data), p=temp_sim_data)
        if math.isnan(LL_temp):
            print('nan LL: no reported infections in this simulation year')
            #print temp_sim_data
        else:
            LL += LL_temp
        #print LL_temp

    return LL

def dirichlet_dengue(raw_data, sim_data):
    num_years = len(raw_data)
    LL = 0.

    for x in range(num_years):
        LL_temp = 0.
        num_strains = len(raw_data[x])
        #print len(raw_data[x])
        raw_nobs = sum(raw_data[x])
        #print type(sim_data[x][1:num_strains+1])
        sim_nobs = sum(sim_data[x][1:num_strains+1])

        LL_temp += gammaln(raw_nobs + 1)
        LL_temp += gammaln(sim_nobs + num_strains)
        LL_temp -= gammaln(raw_nobs + sim_nobs + num_strains)
        for strain in range(num_strains):
            #print strain
            #K.J.S.: Notice that whenever a value is taken from sim_data, the array index is increased by one.  This is because the first
            #	element in the list is supposed to refer to incidence (it may not actually be doing so; all the values are zero)
            LL_temp += gammaln(raw_data[x][strain] + sim_data[x][strain + 1] + 1)  # first column is incidence of all strains
            LL_temp -= gammaln(sim_data[x][strain + 1] + 1)
            LL_temp -= gammaln(raw_data[x][strain] + 1)

        LL_temp /= num_strains
        LL += LL_temp

    return LL

def euclidean_distance(raw_data, sim_data):
    #Convert 2-dimensional array to single dimension
    raw_data = sum(raw_data,[])
    sim_data = sum(sim_data, [])
    num_obs = len(raw_data)
    return math.sqrt(sum([(raw_data[x] - sim_data[x])**2 for x in range(num_obs)]))*-1


def weighted_squares(raw_data, sim_data):
    sim_data = sum(sim_data, [])
    raw_data = sum(raw_data, [])
    num_obs = len(raw_data)
    return math.sqrt(sum([(raw_data[x] - sim_data[x])**2/raw_data[x] for x in range(num_obs)]))*-1
