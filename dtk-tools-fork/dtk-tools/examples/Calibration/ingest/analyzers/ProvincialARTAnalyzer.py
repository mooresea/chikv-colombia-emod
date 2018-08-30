import logging
import pandas as pd
import os, sys
import numpy as np

from dtk.utils.observations.DataFrameWrapper import DataFrameWrapper
from dtk.utils.observations.PopulationObs import PopulationObs
from analyzers.PostProcessAnalyzer import PostProcessAnalyzer

# Plotting
import matplotlib.pyplot as plt
from matplotlib import collections as mc

logger = logging.getLogger(__name__)

class ProvincialARTAnalyzer(PostProcessAnalyzer):

    reference_key = 'On_ART'
    sim_reference_key = 'Sim_On_ART'
    log_float_tiny = np.log( np.finfo(float).tiny )

    def __init__(self, site, weight, **kwargs):
        super(ProvincialARTAnalyzer, self).__init__(**kwargs)

        self.filenames += [ os.path.join('output', 'post_process', 'Number_On_ART.csv') ]

        self.name = self.__class__.__name__
        self.weight = weight
        self.setup = {}

        self.site = site
        self.reference = self.site.reference_data
        self.reference = self.reference.filter(keep_only=[self.reference_key])

    def filter(self, sim_metadata):
        ret = super(ProvincialARTAnalyzer, self).filter(sim_metadata)
        return ret

    def apply(self, parser):
        ret = super(ProvincialARTAnalyzer, self).apply(parser)
        pop_scaling = ret['Pop_Scaling']

        # Read Number on ART:
        sim = parser.raw_data[self.filenames[-1]].copy()

        sim['AgeBin'] = sim['AgeBin'].astype('category', categories=sorted( list(set(sim['AgeBin'])), key=PostProcessAnalyzer.lower_value), ordered=False)
        sim = sim.rename(columns={'Result':self.sim_reference_key, 'Node':'Province'})
        sim[self.sim_reference_key] *= pop_scaling
        sim = sim.set_index('Province')
        sim.rename(self.node_map, inplace=True)
        sim = sim.reset_index()

        stratifiers = ['Year', 'Province', 'Gender', 'AgeBin']
        sim_dfw = DataFrameWrapper(dataframe=sim, stratifiers=stratifiers)
        merged = self.reference.merge(sim_dfw, index=stratifiers, keep_only=[self.reference_key, self.sim_reference_key])

        merged_years = merged.get_years()
        reference_years = self.reference.get_years()

        if reference_years != merged_years:
            raise Exception("[%s] Failed to find all data years (%s) in simulation output (%s)." % (self.name, reference_years, merged_years))

        # If analyzing simulation not generated by itertool, __sample_index__ will not be in tags
        # Instead, use Run_Number
        sample = parser.sim_data.get('__sample_index__')
        if sample is None:
            sample = parser.sim_data.get('Run_Number')

        merged = merged._dataframe
        merged.index.name = 'Index'

        shelve_data = {
            'Data': merged,
            'Sim_Id': parser.sim_id,
            'Sample': sample
        }
        self.shelve_apply( parser.sim_id, shelve_data)

        if self.debug:
            print("size (MB):", sys.getsizeof(shelve_data)/8.0/1024.0)

    def compare_year_gender(self, sample):
    # Note: Might be called extra times by pandas on apply for purposes of "optimization"
    # http://stackoverflow.com/questions/21635915/why-does-pandas-apply-calculate-twice
        #
        log_root_2pi = np.multiply(0.5,np.log(np.multiply(2,np.pi)))
        #
        #
        raw_data = sample['On_ART']
        sim_data = sample[self.sim_reference_key]
        # print('\n---\nraw:\n%s\nsim:\n%s' % (raw_data, sim_data))
        raw_data_variance = np.multiply(sample['On_ART'], 4000*4000) # 400000000)

        # return np.subtract(raw_data,sim_data)
        log_of_gaussian = - log_root_2pi - np.multiply(0.5,np.log(raw_data_variance)) - np.multiply(0.5,np.square(sim_data-raw_data))

        # Scaling

        largest_possible_log_of_gaussian = 0
        largest_possible_log_of_gaussian = largest_possible_log_of_gaussian + (np.multiply(-1, log_root_2pi) - np.multiply(0.5, raw_data_variance) - np.divide(np.multiply(0.5, ((sim_data - raw_data)**2)), raw_data_variance))

        scale_max = 15

        gaussian_ratio = np.divide(scale_max, largest_possible_log_of_gaussian)

        log_of_gaussian = np.multiply(log_of_gaussian, gaussian_ratio)

        # log_of_gaussian = max(log_of_gaussian, self.log_float_tiny)

        return log_of_gaussian

    def compare(self, sample):
        LL = sample.reset_index().groupby(['Year', 'Province', 'Gender']).apply(self.compare_year_gender)
        return (sum(LL.values)*self.weight)

    def combine(self, parsers):
        shelved_data = super(ProvincialARTAnalyzer, self).combine(parsers)

        if shelved_data is not None:
            if self.verbose:
                print('Combine from cache')
            self.data = shelved_data['Data']
            return

        selected = [ self.shelve[str(sim_id)]['Data'] for sim_id in self.sim_ids ]
        keys = [ (self.shelve[str(sim_id)]['Sample'], self.shelve[str(sim_id)]['Sim_Id'])
            for sim_id in self.sim_ids ]

        self.data = pd.concat( selected, axis=0,
                            keys=keys,
                            names=['Sample', 'Sim_Id'] )

        self.data.reset_index(level='Index', drop=True, inplace=True)

        try:
            self.shelve_combine({'Data':self.data})
        except:
            print("shelve_combine didn't work, sorry")

    def cache(self):
        pass

    def uid(self):
        print('UID')
        ''' A unique identifier of site-name and analyzer-name. '''
        return '_'.join([self.site.name, self.name])

    def make_collection(self, d):
        return zip(d['Year'], d['Sim_Prevalence'])

    def plot_agebin(self, data):
        agebin = data['AgeBin'].values[0] # Could be more efficient
        data.set_index('Gender', inplace=True)
        genders = data.index.unique().values.tolist()

        ref = self.reference.reset_index().set_index(['AgeBin', 'Gender']).loc[agebin]

        fig, ax = plt.subplots(1, len(genders), figsize=(16,10), sharey='row', sharex='row')
        for gender, a in zip(genders, ax):
            data_g = data.loc[[gender]]
            ref_g = ref.loc[gender]

            # Color by result?
            data_g_by_sim_id = data_g.groupby('Sim_Id')
            lc = mc.LineCollection( data_g_by_sim_id.apply(self.make_collection), linewidths=0.1, cmap=plt.cm.jet )
            #lc.set_array( data_g_by_sim_id['Results'].apply(lambda z:z)) # <-- Hopefully same order?

            a.add_collection(lc)

            # Use Count to make a poisson confidence interval (noraml approx)
            a.plot(ref_g['Year'], ref_g['NationalPrevalence'], 'k.', ms=25)

            a.autoscale()
            a.margins(0.1)
            a.set_title('%s: %s'%(gender,agebin))
        return fig, agebin


    def finalize(self):
        fn = 'Results_%s.csv'%self.__class__.__name__
        out_dir = os.path.join(self.working_dir, self.basedir, self.exp_id)
        print('--> Writing %s to %s'%(fn, out_dir))
        ProvincialARTAnalyzer.mkdir_p(out_dir)
        results_filename = os.path.join(out_dir, fn)
        self.data.to_csv(results_filename)

        # Call 'compare' ... TODO: Check!
        self.result = self.data.reset_index().groupby(['Sample']).apply(self.compare)

        # Close the shelve file, among other little things.  Can take a long time:
        super(ProvincialARTAnalyzer, self).finalize()
