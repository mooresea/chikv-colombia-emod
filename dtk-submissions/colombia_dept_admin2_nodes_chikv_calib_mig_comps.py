# Execute directly: 'python example_calibration.py'
# or via the calibtool.py script: 'calibtool run example_calibration.py'
# Execute directly: 'python example_calibration.py'
# or via the calibtool.py script: 'calibtool run example_calibration.py'
import math
import os
import json
import sys

from calibtool.CalibManager import CalibManager
from calibtool.algorithms.OptimTool import OptimTool
from calibtool.plotters.LikelihoodPlotter import LikelihoodPlotter
from scipy.special import gammaln

#sys.path.append(os.path.abspath('Single_Node_Sites'))
from outbreakindividualdengue import add_OutbreakIndivisualDengue
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.study_sites import configure_site
from scipy.stats import uniform
from simtools.SetupParser import SetupParser
from ColombiaChikvTSSite_Department_admin2_nodes import ColombiaTSSite
import numpy as np
import scipy as sp

site_name = sys.argv[1]

SetupParser.default_block = "HPC"
with SetupParser.TemporarySetup("LOCAL") as setp:
    input_dir = setp.get('input_root')

# @@ should be imported from elsewhere. PENDING until decided on structure for this.
def foi_time(times, widths, heights, days_sim, node_cnt):

    foi_out = [np.repeat(0, repeats=len(days_sim)) for ii in range(0,node_cnt)]

    for ss in range(0, node_cnt):
        for yy in range(0, len(times[ss])):
            times[ss][yy] = times[ss][yy]*len(days_sim)
            foi_out[ss] = \
                foi_out[ss] + \
                heights[ss][yy] * \
                sp.stats.norm(times[ss][yy], widths[ss][yy]).pdf(days_sim)

    return foi_out

cb = DTKConfigBuilder.from_defaults('DENGUE_SIM')  # should be from_files (check what cls means)
#configure_site(cb, 'Puerto_Rico')
cb.set_param('Koppen_Filename',"")
climate_file_path = 'D:\\Eradication/InputDataFiles/Colombia_dept_nodes_admin2/%s/' % site_name
cb.set_input_files_root(climate_file_path)
cb.set_experiment_executable('D:\\Eradication\\Executable\\Eradication.exe')
cb.set_dll_root('D:\\Eradication\\bin')
cb.set_param('Air_Temperature_Filename','%s_admin2_air_temperature_daily.bin' % site_name)
cb.set_param('Land_Temperature_Filename','%s_admin2_air_temperature_daily.bin' % site_name)
cb.set_param('Rainfall_Filename','%s_admin2_rainfall_daily.bin' % site_name)
cb.set_param('Relative_Humidity_Filename','%s_admin2_humidity_daily.bin' % site_name)
demo_file_name = "%s_admin2_nodes_demographics.json" % site_name
mig_file_name = "%s_admin2_nodes_regional_migration_hi.bin" % site_name
cb.set_param('Demographics_Filenames',[demo_file_name])
cb.set_param('Regional_Migration_Filename',mig_file_name)


## Get node ids
demographics_path = os.path.join(climate_file_path, cb.get_param('Demographics_Filenames')[0])
fp = open(demographics_path,'rb')
demog = json.load(fp)
node_cnt = demog['Metadata']['NodeCount']
node_ids = []
node_pops = []
for i in range(0,len(demog['Nodes'])):
    node_ids.append(demog['Nodes'][i]['NodeID'])
    node_pops.append(demog['Nodes'][i]['NodeAttributes']['InitialPopulation'])
#Create dictionary to store node_ids and node_pop
#node_dict=dict(zip(node_ids,node_pop))

cb.enable('Spatial_Output')
cb.set_param('Spatial_Output_Channels',["New_Reported_Infections","Population"])
#set_larval_habitat(cb,{"aegypti":{'TEMPORARY_RAINFALL':pow(10,13)}})
cb.config["parameters"]["Strain1_Background_Exposure"]=0.00
#cb.config["parameters"]["Strain2_Background_Exposure"]=0.001
#cb.config["parameters"]["Strain4_Background_Exposure"]=0.001
scale_factor = 1.0
max_pop = 20000
cb.set_param("Base_Population_Scale_Factor",scale_factor)
cb.set_param('x_Temporary_Larval_Habitat', scale_factor)
cb.params['Vector_Sampling_Type'] = 'VECTOR_COMPARTMENTS_NUMBER'
cb.config["parameters"]["Max_Node_Population_Samples"]=max_pop
cb.set_param("Individual_Sampling_Type", "ADAPTED_SAMPLING_BY_POPULATION_SIZE")
scaled_pop = np.asarray(node_pops) * scale_factor
nodes_pop_scaled=dict(zip(node_ids,scaled_pop))
cb.set_param('Enable_Demographics_Reporting',0)

##Set migration
cb.set_param('Migration_Model','FIXED_RATE_MIGRATION')
cb.set_param('Migration_Pattern','SINGLE_ROUND_TRIPS')
cb.set_param('Enable_Regional_Migration',1)
cb.set_param('Regional_Migration_Roundtrip_Duration',1.0)
cb.set_param('Regional_Migration_Roundtrip_Probability', 1)
cb.set_param('x_Regional_Migration',1.0)

##Get sampling pop size for use with importations
for node_id in range(node_cnt):
    if scaled_pop[node_id] > max_pop:
        scaled_pop[node_id] = max_pop
nodes_scaled=dict(zip(node_ids,scaled_pop))
####
#### CHIKV specific parameter values
####
##Infectiousness
cb.config["parameters"]["Vector_Species_Params"]["aegypti"]["Infected_Arrhenius_1"] = 9.47*pow(10,12)
cb.set_param('Infectiousness_Asymptomatic_Naive_1',0.547)
cb.set_param('Infectiousness_Asymptomatic_Naive_2',3.2564)
cb.set_param('Infectiousness_Asymptomatic_Naive_3',1.488862)
cb.set_param('Infectiousness_Asymptomatic_Reinfected_1',0.547)
cb.set_param('Infectiousness_Asymptomatic_Reinfected_2',3.2564)
cb.set_param('Infectiousness_Asymptomatic_Reinfected_3',1.488862)
cb.set_param('Infectiousness_Symptomatic_Naive_1',0.547)
cb.set_param('Infectiousness_Symptomatic_Naive_2',3.2564)
cb.set_param('Infectiousness_Symptomatic_Naive_3',1.488862)
cb.set_param('Infectiousness_Symptomatic_Reinfected_1',0.547)
cb.set_param('Infectiousness_Symptomatic_Reinfected_2',3.2564)
cb.set_param('Infectiousness_Symptomatic_Reinfected_3',1.488862)

##Symptomatic Rates
symptom_rate = 0.72
report_rate = 0.08 * symptom_rate
cb.set_param('Symptomatic_Probability_Naive',symptom_rate)
cb.set_param('Symptomatic_Probability_Reinfected',symptom_rate)
cb.set_param('Symptomatic_Probability_Postsecondary',symptom_rate)
cb.set_param('Symptomatic_Reporting_Rate',report_rate)

##Incubation period and reporting delays
cb.set_param('Incubation_Period_Log_Mean',log(3))
cb.set_param('Incubation_Period_Log_Width',log(1.04))
cb.set_param('Reporting_Period_Log_Mean',log(7))
cb.set_param('Reporting_Period_Log_Width',log(2))
###Simulation parameters
sim_length = 159*7 #Length of Colombia timeseries
cb.set_param('Simulation_Duration',sim_length) #86*7)
cb.campaign["Campaign_Name"] = "Campaign - Outbreak"
# add_OutbreakIndivisualDengue(cb, 520, {},0.0005, 'Strain_1', [])
# add_OutbreakIndivisualDengue(cb, 580, {},0.0005, 'Strain_1', [])
# add_OutbreakIndivisualDengue(cb, 640, {},0.0005, 'Strain_1', [])
cb.params['.logLevel_JsonConfigurable'] = 'WARNING'

sites = [ColombiaTSSite(site_name)]
plotters = [LikelihoodPlotter(True)]

def create_param(name, guess, min, max):
    return { 'Name': name, 'Dynamic': True, 'Guess': guess, 'Min': min, 'Max': max }

params = []

# Add TEMPORARY_RAINFALL
params.append(create_param('TEMPORARY_RAINFALL', 10, 0.1, 100))
params.append(create_param('Temporary_Habitat_Decay_Factor', 0.05, 0.005, 0.1))
params.append(create_param('x_Regional_Migration', 0.1, 0.01, 1.0))

#Add the Import_Times, Import_Heights, and Import_Widths by node and time
#For now only allow one introduction per node
#Code written so multiple introductions can be simulated
intro_cnt=1
###Set number of nodes that importation will occur in
import_node_cnt=1
##Now select which node(s) for importation
days_import = range(1,sim_length)
for i in range(0, import_node_cnt):
    ##Select node based on pop size
    i_node = sorted(nodes_pop_scaled,key=nodes_pop_scaled.__getitem__,reverse=True)[i]
    for j in range(0, intro_cnt):
        ## First case in Colombia reported in week 74 (day 511 or later)
        ## First case in Caribbean reported in late 2013
        ## Will assume cases in Colombia didn't appear until after first St Martin case (day 365 - 0.33 into timeseries)
        ## Useful since mosquito population needs 1 year to rampup and earlier intros could give misleading results
        params.append(create_param('Import_Times_%s_%s' % (i_node,j), 0.5,0.33,1))
        params.append(create_param('Import_Heights_%s_%s' % (i_node,j),10,0.001,100))
        params.append(create_param('Import_Widths_%s_%s' % (i_node,j),5,1,50))

def sample_point_fn(cb, params_dict):
    """
    A simple example function that takes a list of sample-point values
    and sets parameters accordingly using the parameter names from the prior.
    Note that more complicated logic, e.g. setting campaign event coverage or habitat abundance by species,
    can be encoded in a similar fashion using custom functions rather than the generic "set_param".
    """
    cb.config["parameters"]["Vector_Species_Params"]["aegypti"]["Larval_Habitat_Types"]["TEMPORARY_RAINFALL"] = params_dict['TEMPORARY_RAINFALL']*pow(10,10)
    cb.config["parameters"]['Temporary_Habitat_Decay_Factor'] = params_dict['Temporary_Habitat_Decay_Factor']
    cb.config["parameters"]['x_Regional_Migration'] = params_dict['x_Regional_Migration']

    par_times_import = []
    par_widths_import = []
    par_heights_import = []
    for i in range(0, import_node_cnt):
        i_node = sorted(nodes_pop_scaled,key=nodes_pop_scaled.__getitem__,reverse=True)[i]
        par_times_import.append([(params_dict['Import_Times_%s_%s' % (i_node,str(yy))]) for yy in range(0, intro_cnt)])
        par_widths_import.append([(params_dict['Import_Widths_%s_%s' % (i_node,str(yy))]) for yy in range(0, intro_cnt)])
        par_heights_import.append([(params_dict['Import_Heights_%s_%s' % (i_node,str(yy))]) for yy in range(0, intro_cnt)])

    # generate time-varying force of infection patterns by serotype
    foi_times = foi_time(par_times_import, par_widths_import, par_heights_import, days_import, import_node_cnt)

    # simulate importations from the force of importation patterns
    importations = []

    for i in range(0, len(foi_times)):
        importations.append(np.random.poisson(foi_times[i]))

    #Importations are a proportion of (scaled) population = number of importations[node][time] / population_size[node]
    for tt in range(0, sim_length - 1):
        for node_id in range(import_node_cnt):
            i_node = sorted(nodes_pop_scaled,key=nodes_pop_scaled.__getitem__,reverse=True)[i]
            if importations[node_id][tt] > 0.0:
                add_OutbreakIndivisualDengue(cb, tt, {}, importations[node_id][tt] / nodes_scaled[i_node], 'Strain_1',
                                             [i_node])

    return params_dict

# Just for fun, let the numerical derivative baseline scale with the number of dimensions
volume_fraction = 0.05   # desired fraction of N-sphere area to unit cube area for numerical derivative (automatic radius scaling with N)
num_params = len( [p for p in params if p['Dynamic']] )
r = math.exp( 1/float(num_params)*( math.log(volume_fraction) + gammaln(num_params/2.+1) - num_params/2.*math.log(math.pi) ) )

optimtool = OptimTool(params, lambda p: p,
    mu_r = r,           # <-- radius for numerical derivatve.  CAREFUL not to go too small with integer parameters
    sigma_r = r/10.,    # <-- stdev of radius
    center_repeats = 1, # <-- Number of times to replicate the center (current guess).  Nice to compare intrinsic to extrinsic noise
    samples_per_iteration = 100 # <-- Samples per iteration, includes center repeats.  Actual number of sims run is this number times number of sites.
)

calib_name = "CHIKV_Calib_Habitat_Imports_admin2nodes_Migration" \
             "_%s" % site_name
calib_manager = CalibManager(name=calib_name,
                             config_builder=cb,
                             map_sample_to_model_input_fn=sample_point_fn,
                             sites=sites,
                             next_point=optimtool,
                             sim_runs_per_param_set=4,
                             max_iterations=10,
                             plotters=plotters)

#run_calib_args = {'selected_block': "LOCAL"}
run_calib_args = {}

if __name__ == "__main__":
    SetupParser.init(selected_block=SetupParser.default_block)
    #run_calib_args.update(dict(location='HPC'))
    calib_manager.cleanup()
    calib_manager.run_calibration() #(**run_calib_args)




