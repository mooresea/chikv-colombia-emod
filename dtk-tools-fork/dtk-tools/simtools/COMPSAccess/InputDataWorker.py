import json


class InputDataWorker:
    def __init__(self, demographics_file_path, wo_output_path, project_info="IDM-Democratic_Republic_of_the_Congo",
                 include_non_pop=True, shape_id="", resolution="30", parameters=('tmean', 'humid', 'rain'),
                 start_year='2009', num_years="4", nan_check=True, idRef='Gridded world grump2.5arcmin', project_root='v2017'):

        self.demographics_file_path = demographics_file_path
        self.wo_output_path = wo_output_path
        self.work_item_type = "InputDataWorker"
        self.node_list = []
        self.project_info = project_info
        self.region = ""
        self.include_non_pop = include_non_pop
        self.shape_id = shape_id
        self.resolution = resolution
        self.parameters = parameters
        self.start_year = start_year
        self.num_years = num_years
        self.nan_check = nan_check
        self.project_root = project_root
        self.migration = False
        self.id_reference = idRef
        self.path = ''

    def wo_2_dict(self):
        wo = {'WorkItem_Type': self.work_item_type,
              'Project': self.project_info,
              'Region': self.region,
              'IncludeNonPop': self.include_non_pop,
              'Resolution': self.resolution,
              'Parameters': self.parameters,
              'StartYear': self.start_year,
              'NumYears': self.num_years,
              'NaNCheck': self.nan_check,
              'Migration': self.migration,
              'IdReference': self.id_reference,
              'Mode': 'upload',
              'ProjectRoot':self.project_root
              }

        # add nodes from demographics file
        with open(self.demographics_file_path, 'r') as demo_f:
            demographics = json.load(demo_f)

            # generate node list for work order
            for node in demographics['Nodes']:
                self.node_list.append(node['NodeID'])

        # wo['NodeList'] = self.node_list

        return wo

    def wo_2_json(self):
        with open(self.wo_output_path, 'w') as wo_f:
            json.dump(self.wo_2_dict(), wo_f, indent=3)
