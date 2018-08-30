import json
import os
from calibtool.CalibSite import CalibSite

from CHIKVTimeSeriesAnalyzer_combinednodes import CHIKVTimeSeriesAnalyzer  # @@this is not imported in PuertoRicoSite.py. Why?

class ColombiaTSSite(CalibSite):
    def __init__(self,site_name):
        # type: () -> object
        super(ColombiaTSSite, self).__init__(site_name)

    def get_reference_data(self, reference_type):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        input_file_name = self.name + '_2013-2016_timeseries.json'
        with open(os.path.join(dir_path,'inputs',input_file_name)) as data_file:
            ref = json.load(data_file)

        return ref[reference_type]

    def get_analyzers(self):
        return [CHIKVTimeSeriesAnalyzer(self)]

    def get_setup_functions(self):
        """
        Derived classes return a list of functions to apply site-specific modifications to the base configuration.
        These are combined into a single function using the SiteFunctions helper class in the CalibSite constructor.
        """
        return []
