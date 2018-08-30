import os
from simtools.SetupParser import SetupParser
from simtools.Utilities.General import init_logging

logger = init_logging('ExperimentManager')


class ExperimentManagerFactory(object):
    @staticmethod
    def _factory(type):
        if type == 'LOCAL':
            from simtools.ExperimentManager.LocalExperimentManager import LocalExperimentManager
            return LocalExperimentManager
        if type == 'HPC':
            from simtools.ExperimentManager.CompsExperimentManager import CompsExperimentManager
            return CompsExperimentManager
        raise Exception("ExperimentManagerFactory location argument should be either 'LOCAL' or 'HPC'.")

    @classmethod
    def from_experiment(cls, experiment, config_builder=None):
        # If a string is passed -> get the experiment from DB
        if isinstance(experiment, str):
            from simtools.Utilities.Experiments import retrieve_experiment
            experiment = retrieve_experiment(experiment)

        logger.debug("Factory - Creating ExperimentManager for experiment %s pid: %d location: %s" % (experiment.id, os.getpid(), experiment.location))
        manager_class = cls._factory(type=experiment.location)
        manager = manager_class(experiment=experiment, config_builder=config_builder)

        return manager

    @classmethod
    def from_model(cls, model_file, location='LOCAL', config_builder=None):
        logger.debug('Factory - Initializing %s ExperimentManager from: %s', location, model_file)
        return cls._factory(location)(experiment=None, config_builder=config_builder)

    @classmethod
    def from_setup(cls):
        location = SetupParser.get('type')
        logger.debug('Factory - Initializing %s ExperimentManager from parsed setup' % location)
        logger.warning('ExperimentManagerFactory.from_setup is deprecated. Use ExperimentManagerFactory.init() instead or '
                       'ExperimentManagerFactory.from_cb(config_builder) if you have a config_builder')
        return cls._factory(location)(experiment=None, config_builder=None)

    @classmethod
    def init(cls):
        location = SetupParser.get('type')
        logger.debug('Factory - Initializing %s ExperimentManager' % location)
        return cls._factory(location)(experiment=None, config_builder=None)

    @classmethod
    def from_cb(cls, config_builder=None):
        location = SetupParser.get('type')
        logger.debug('Factory - Initializing %s ExperimentManager from config_builder' % location)
        return cls._factory(location)(experiment=None, config_builder=config_builder)


