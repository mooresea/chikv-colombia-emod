import logging
import os

import pandas as pd
import LL_calculators_dengue
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from calibtool import LL_calculators
from calibtool.analyzers.BaseComparisonAnalyzer import BaseComparisonAnalyzer
from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer

logger = logging.getLogger(__name__)


class CHIKVTimeSeriesAnalyzer(BaseComparisonAnalyzer):
    filenames = ['output/InsetChart.json','config.json']
    config_filename =['config.json']

    x = 'Week'
    y = 'Reported Cases'

    data_group_names = ['sample', 'sim_id', 'channel']

    def __init__(self,site, weight=1,):
        super(CHIKVTimeSeriesAnalyzer, self).__init__(site, weight)
        self.name = "IncidenceAnalyze"
        self.reference = site.get_reference_data('weekly_reported_cases')

    def set_site(self, site):
        '''
        Get the reference data that this analyzer needs from the specified site.
        '''

        self.site = site
        self.reference = self.site.reference_data["weekly_reported_cases"]

    def filter(self, sim_metadata):
        '''
        This analyzer only needs to analyze simulations for the site it is linked to.
        N.B. another instance of the same analyzer may exist with a different site
             and correspondingly different reference data.
        '''
        return sim_metadata.get('__site__', False) == self.site.name

    def apply(self, parser):
        '''
        Extract data from output data and accumulate in same bins as reference.
        '''
        #First file is infections by node
        data = parser.raw_data[self.filenames[0]]
        inf_data=data['Channels']['New Reported Infections']['Data']
        infection_data = pd.Series(data['Channels']['New Reported Infections']['Data'])
        #Get node IDs for matching with reference data
        #node_ids = parser.raw_data[self.filenames[0]]['nodeids']

        #Set node names to distinguish b/w population and infections when datasets combined into single dataframe
        # inf_node_names = []
        # pop_node_names = []
        #
        # for i in range(0,len(node_ids)):
        #     inf_node_names.append('%s_infections' % node_ids[i])
        #     pop_node_names.append('%s_populations' % node_ids[i])

        population = pd.Series(data['Channels']['Statistical Population']['Data'])

        # incidence = []
        # for i in range( len(new_infections_weekly)):
        #     #curr = new_infections_weekly[i] / population_weekly[i]
        #     curr = new_infections_weekly[i] * 100.0 ## scale factor between sample_population size and Department pop size
        #     incidence.append(curr)

       # print incidence

        #Convert to dataframes for combining samples/simulations
        # infection_df = pd.DataFrame(data=infection_data,columns=inf_node_names,index=range(1,len(infection_data)+1))
        # infection_df.index.name='day'
        # population_df = pd.DataFrame(data=population,columns=pop_node_names,index=range(1,len(infection_data)+1))
        # population_df.index.name='day'

        ###Aggregate daily simulation data to weekly

        infection_weekly = infection_data.groupby(infection_data.index / 7).sum()
        population_weekly = population.groupby(population.index / 7).mean()

        channel_data = pd.DataFrame({'Infections': infection_weekly,
                                    'Population': population_weekly},
                                    index=range(0,len(infection_weekly)))
        channel_data.index.name = 'week'
        #channel_data = infection_weekly.merge(population_weekly, left_index=True,right_index=True)

        ##Retrieve simulation population scaling factor so that infections and population match full population
        scale_factor = 1.0 / parser.raw_data[self.filenames[1]]['parameters']['Base_Population_Scale_Factor']
        channel_data = channel_data.multiply(scale_factor)
        channel_data.sample = parser.sim_data.get('__sample_index__')
        channel_data.sim_id = parser.sim_id

        return channel_data

    def combine(self, parsers):
        '''
        Combine the simulation data into a single table for all analyzed simulations.
        '''

        selected = [p.selected_data[id(self)] for p in parsers.values() if id(self) in p.selected_data]
        combined = pd.concat(selected, axis=1,
                             keys=[(d.sample, d.sim_id) for d in selected],
                             names=self.data_group_names)
        stacked = combined.stack(['sample','sim_id'])
        self.data = stacked.groupby(level=['sample', 'week']).mean()
       # self.data = combined.groupby(level=['sample', 'channel'], axis=1).mean()

    def compare(self, sample):
        '''
        Assess the result per sample, in this case the likelihood
        comparison between simulation and reference data.
        '''
        ##Convert observed data to DataFrame so that same manipulations can be performed on observed and simulated data
        # if(type(self.reference['Node IDs']) is int):
        #     rep_data_df = pd.DataFrame(data=self.reference["Reported Cases"],
        #                                columns=[str(self.reference['Node IDs'])])
        # else:
        #     rep_data_df = pd.DataFrame(data=self.reference["Reported Cases"],
        #                                columns=map(str, self.reference['Node IDs']))
        #Make sure data sorted by node ID (as string)
        #rep_data_df = rep_data_df[sorted(rep_data_df.columns)]
        rep_data=self.reference["Reported Cases"]

        #Ignore time before first possible introduction

        ##Extract population and infections into separate data frames
        pop_sample = sample['Population']
        #print pop_sample.values.tolist()
        inf_sample = sample['Infections']
        #print inf_sample.values.tolist()
        #print rep_data

        self.subplot.plot(inf_sample.values.tolist())
        return LL_calculators_dengue.betaBinomial_dengue_simple(rep_data,
                                                inf_sample.values.tolist(),pop_sample.values.tolist())
        #return LL_calculators_dengue.euclidean_distance(self.reference["Reported Cases"],
        #                                         sample.values.tolist())

    def finalize(self):
        '''
        Calculate the output result for each sample.
        '''
        self.figure = Figure()
        self.subplot = self.figure.add_subplot(111)
        self.subplot.plot(self.reference["Reported Cases"])

        self.result = self.data.groupby(level='sample').apply(self.compare)
        #       self.result = self.data.apply(self.compare)
        logger.debug(self.result)

        canvas = FigureCanvasAgg(self.figure)
        canvas.print_figure(os.path.join(self.working_dir, 'fig.png'), dpi=80)

    def cache(self):
        '''
        Return a cache of the minimal data required for plotting sample comparisons
        to reference comparisons.
        '''

        return {}

    def uid(self):
        ''' A unique identifier of site-name and analyzer-name. '''
        return '_'.join([self.site.name, self.name])

    def plot_comparison(cls, fig, data, **kwargs):
        pass
