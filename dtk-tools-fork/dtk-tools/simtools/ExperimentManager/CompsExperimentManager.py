import os
from multiprocessing import Process

from COMPS.Data import Experiment, Configuration, Priority, Suite

from simtools.ExperimentManager.BaseExperimentManager import BaseExperimentManager
from simtools.SetupParser import SetupParser
from simtools.SimulationCreator.COMPSSimulationCreator import COMPSSimulationCreator
from simtools.Utilities.COMPSCache import COMPSCache
from simtools.Utilities.COMPSUtilities import experiment_is_running, COMPS_login, get_semaphore, \
    get_simulation_by_id
from simtools.Utilities.General import init_logging, timestamp

logger = init_logging("COMPSExperimentManager")


class CompsExperimentManager(BaseExperimentManager):
    """
    Extends the LocalExperimentManager to manage DTK simulations through COMPSAccess wrappers
    e.g. creation of Simulation, Experiment, Suite objects
    """
    location = 'HPC'

    def __init__(self, experiment=None, config_builder=None):
        # Ensure we use the SetupParser environment of the experiment if it already exists
        super().__init__(experiment, config_builder)
        temp_block = experiment.selected_block if experiment else SetupParser.selected_block
        temp_path = experiment.setup_overlay_file if experiment else None
        self.endpoint = experiment.endpoint if experiment else None

        with SetupParser.TemporarySetup(temporary_block=temp_block, temporary_path=temp_path) as setup:
            self.comps_sims_to_batch = int(setup.get(parameter='sims_per_thread'))
            self.endpoint = self.endpoint or setup.get(parameter='server_endpoint')
            COMPS_login(self.endpoint)

        self.asset_service = True
        self.runner_thread = None
        self.sims_to_create = []
        self.commissioners = []
        self.save_semaphore = get_semaphore()

        # If we pass an experiment, retrieve it from COMPS
        if self.experiment:
            self.comps_experiment = COMPSCache.experiment(self.experiment.exp_id)

    def get_simulation_creator(self, work_queue):
        return COMPSSimulationCreator(config_builder=self.config_builder,
                                      initial_tags=self.exp_builder.tags,
                                      work_queue=work_queue,
                                      experiment=self.experiment,
                                      cache=self.cache,
                                      save_semaphore=self.save_semaphore,
                                      comps_experiment=self.comps_experiment)

    @staticmethod
    def create_suite(suite_name):
        suite = Suite(suite_name)
        suite.save()

        return str(suite.id)

    MAX_SUBDIRECTORY_LENGTH = 50 - len(timestamp()) - 1 # avoid maxpath issues on COMPS

    def create_experiment(self, experiment_name, experiment_id=None, suite_id=None):
        # Also create the experiment in COMPS to get the ID
        COMPS_login(SetupParser.get('server_endpoint'))
        experiment_name = self.clean_experiment_name(experiment_name)
        subdirectory = experiment_name[0:self.MAX_SUBDIRECTORY_LENGTH] + '_' + timestamp()
        config = Configuration(
            environment_name=SetupParser.get('environment'),
            simulation_input_args=self.commandline.Options,
            working_directory_root=os.path.join(SetupParser.get('sim_root'), subdirectory),
            executable_path=self.commandline.Executable,
            node_group_name=SetupParser.get('node_group'),
            maximum_number_of_retries=int(SetupParser.get('num_retries')),
            priority=Priority[SetupParser.get('priority')],
            min_cores=self.config_builder.get_param('Num_Cores', 1) if self.config_builder else 1,
            max_cores=self.config_builder.get_param('Num_Cores', 1) if self.config_builder else 1,
            exclusive=self.config_builder.get_param('Exclusive', False) if self.config_builder else False
        )

        e = Experiment(name=experiment_name,
                       configuration=config,
                       suite_id=suite_id)

        # Add tags if present
        if self.experiment_tags: e.set_tags(self.experiment_tags)

        e.save()

        # Store in our object
        self.comps_experiment = e

        # Also add it to the cache
        COMPSCache.add_experiment_to_cache(e)

        # Create experiment in the base class
        super(CompsExperimentManager, self).create_experiment(experiment_name,  str(e.id), suite_id)

        # Set some extra stuff
        self.experiment.endpoint = self.endpoint

    def create_simulation(self):
        files = self.config_builder.dump_files_to_string()

        # Create the tags and append the environment to the tag
        tags = self.exp_builder.metadata
        tags['environment'] = SetupParser.get('environment')
        tags.update(self.exp_builder.tags if hasattr(self.exp_builder, 'tags') else {})

        # Add the simulation to the batch
        self.sims_to_create.append({'name': self.config_builder.get_param('Config_Name'), 'files':files, 'tags':tags})

    def commission_simulations(self):
        """
        Launches an experiment and its associated simulations in COMPS
        :param states: a multiprocessing.Queue() object for simulations to use for updating their status
        :return: The number of simulations commissioned.
        """
        from simtools.SimulationRunner.COMPSRunner import COMPSSimulationRunner
        if self.experiment and (not self.runner_thread or not self.runner_thread.is_alive()):
            logger.debug("Commissioning simulations for COMPS experiment: %s" % self.experiment.id)
            self.runner_thread = Process(target=COMPSSimulationRunner, args=(self.experiment, self.comps_experiment))
            self.runner_thread.daemon = True
            self.runner_thread.start()
            return len(self.experiment.simulations)
        else:
            return 0

    def cancel_experiment(self):
        super(CompsExperimentManager, self).cancel_experiment()
        COMPS_login(self.endpoint)
        if self.comps_experiment and experiment_is_running(self.comps_experiment):
            self.comps_experiment.cancel()

    def hard_delete(self):
        """
        Delete data for experiment and marks the server entity for deletion.
        """
        # Mark experiment for deletion in COMPS.
        COMPS_login(self.endpoint)
        self.comps_experiment.delete()

        # Delete in the DB
        try:
            from simtools.DataAccess.DataStore import DataStore
            DataStore.delete_experiment(self.experiment)
        except:
            pass

    def kill_simulation(self, simulation):
        s = get_simulation_by_id(simulation)
        s.cancel()

    def merge_tags(self, additional_tags):
        if self.comps_experiment:
            self.comps_experiment.merge_tags(additional_tags)