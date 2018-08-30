def in_common_run_and_catalyst_arguments(subparser):
    subparser.add_argument(dest='config_name', default=None,
                           help='Name of configuration python script for custom running of simulation.')
    subparser.add_argument('--priority', default=None, help='Specify priority of COMPS simulation (only for HPC).')
    subparser.add_argument('--node_group', default=None, help='Specify node group of COMPS simulation (only for HPC).')
    subparser.add_argument('-q', '--quiet', action='store_true', help='Runs quietly.')
    return subparser

# 'dtk run' options
def populate_run_arguments(subparsers, func):
    parser_run = subparsers.add_parser('run', help='Run one or more simulations configured by run-options.')
    parser_run = in_common_run_and_catalyst_arguments(parser_run)
    parser_run.add_argument('-b', '--blocking', action='store_true', help='Block the thread until the simulations are done.')
    parser_run.set_defaults(func=func)


# 'dtk ini' options
def populate_ini_arguments(subparsers, func):
    parser_ini = subparsers.add_parser('ini', help='Check  simtools.ini file.')
    parser_ini.add_argument('-o', '--open', action='store_true', help='Open and display file.')
    parser_ini.set_defaults(func=func)


# 'dtk run' options for catalyst, after all, it runs and then does a set-piece analysis
def populate_catalyst_arguments(subparsers, func):
    import catalyst_report
    import os

    parser_catalyst = subparsers.add_parser('catalyst', help='Run a timestep or population-scaling sweep to evaluate '
                                                        'model performance.')
    parser_catalyst = in_common_run_and_catalyst_arguments(parser_catalyst)
    parser_catalyst.add_argument('-s', '--sweep_type', dest='sweep_type',
                                 choices=['timestep', 'popsampling'], default='popsampling',
                                 help='The type of performance sweep to run and report on (Default: popsampling).')
    parser_catalyst.add_argument('-m', '--sweep_method', dest='sweep_method', type=str, default=None,
                                 help='The sweeping method to use (Default: depends on sweep_type).')
    parser_catalyst.add_argument('-r', '--report', dest='report_type', type=str, default=None,
                                 help='The type of report to generate '
                                      '(Default: determined by Simulation_Type in config.json).')
    parser_catalyst.add_argument('-id', '--id', dest='experiment_id', default=None,
                                help='Experiment ID to generate a report for. No new simulations are run if provided '
                                     '(Default: run simulations and generate a report).')
    parser_catalyst.add_argument('--start_step', default=None, type=int, dest='step_from',
                                 help="Starting time step for analysis (Default: near beginning of simulations).")
    parser_catalyst.add_argument('--end_step', default=None, type=int, dest='step_to',
                                 help="Ending time step for analysis (Default: end of simulations).")
    parser_catalyst.add_argument('--raw_data', default=False, action='store_true',
                                 help="Saves raw simulation data into a raw_data.csv file. This option may noticeably "
                                      "increase the report generation time (Default: False).")
    parser_catalyst.add_argument('-l', '--report_label', default=None, type=str,
                                 help='Additional descriptive label to attach to reporting directory (Default: None)')


    parser_catalyst.add_argument('--sweep_definitions', default=None, type=str,
                                 help='A JSON file that defines the available sweeps for use (Default: built-in JSON '
                                      'at %s' % os.path.join(catalyst_report.__path__[0], '{pop_sampling.json, time_steps.json} )'))
    parser_catalyst.add_argument('--report_definitions', default=None, type=str,
                                 help='A JSON file that defines the available report types '
                                      '(Default: %s)' % os.path.join(catalyst_report.__path__[0], 'reports.json'))

    parser_catalyst.set_defaults(func=func)


def populate_diskspace_arguments(subparsers, func):
    parser_status = subparsers.add_parser('diskspace', help='Report diskspace of for users.')
    parser_status.add_argument('-u', '--users', nargs='*', help='List of users.', dest='users')
    parser_status.add_argument('-t', '--top', nargs='?', type=int, help='Get given number recent experiment list.', dest='top')
    parser_status.add_argument('-s', '--save', action='store_true', help='Save results to a file.', dest='save')
    parser_status.add_argument('-r', '--refresh', action='store_true', help='Refresh data from COMPS.', dest='refresh')
    parser_status.set_defaults(func=func)


# 'dtk status' options
def populate_status_arguments(subparsers, func):
    parser_status = subparsers.add_parser('status', help='Report status of simulations in experiment specified by ID or name.')
    parser_status.add_argument(dest='expId', default=None, nargs='?', help='Experiment ID or name.')
    parser_status.add_argument('-r', '--repeat', action='store_true', help='Repeat status check until job is done processing.')
    parser_status.add_argument('-a', '--active', action='store_true', help='Get the status of all active experiments (mutually exclusive to all other options).')
    parser_status.set_defaults(func=func)


# 'dtk list' options
def populate_list_arguments(subparsers, func):
    parser_list = subparsers.add_parser('list', help='Report recent 20 list of simulations in experiment.')
    parser_list.add_argument(dest='exp_name', default=None, nargs='?', help='Experiment name.')
    parser_list.add_argument('-n', '--number', help='Get given number recent experiment list', dest='limit')
    parser_list.set_defaults(func=func)


# 'dtk kill' options
def populate_kill_arguments(subparsers, func):
    parser_kill = subparsers.add_parser('kill', help='Kill most recent running experiment specified by ID or name.')
    parser_kill.add_argument(dest='expId', default=None, nargs='?', help=' Experiment ID or name.')
    parser_kill.add_argument('-id', dest='simIds', default=None, nargs='+', help='Process or job IDs of simulations to kill.')
    parser_kill.set_defaults(func=func)


# 'dtk exterminate' options
def populate_exterminate_arguments(subparsers, func):
    parser_exterminate = subparsers.add_parser('exterminate', help='Kill ALL experiments matched by ID or name.')
    parser_exterminate.add_argument(dest='expId', default=None, nargs='?', help=' Experiment ID or name.')
    parser_exterminate.set_defaults(func=func)


# 'dtk delete' options
def populate_delete_arguments(subparsers, func):
    parser_delete = subparsers.add_parser('delete', help='Delete the most recent experiment specified by ID or name.')
    parser_delete.add_argument(dest='expId', default=None, nargs='?', help=' Experiment ID or name.')
    parser_delete.set_defaults(func=func)


# 'dtk clean' options
def populate_clean_arguments(subparsers, func):
    parser_clean = subparsers.add_parser('clean', help='Hard deletes ALL experiments in {current_dir}\simulations matched by ID or name.')
    parser_clean.add_argument(dest='expId', default=None, nargs='?', help=' Experiment ID or name.')
    parser_clean.set_defaults(func=func)


# 'dtk link' options
def populate_link_arguments(subparsers, func):
    parser_exterminate = subparsers.add_parser('link', help='Open browser to COPMS experiment matched by ID or name.')
    parser_exterminate.add_argument(dest='Id', default=None, nargs='?', help=' Experiment ID (or name) or Simulation ID.')
    parser_exterminate.set_defaults(func=func)


# 'dtk stdout' options
def populate_stdout_arguments(subparsers, func):
    parser_stdout = subparsers.add_parser('stdout', help='Print stdout from first simulation in selected experiment.')
    parser_stdout.add_argument(dest='expId', default=None, nargs='?', help=' Experiment ID or name.')
    parser_stdout.add_argument('-e', '--error', action='store_true', help='Print stderr instead of stdout.')
    parser_stdout.add_argument('--failed', action='store_true', help='Get the stdout for the first failed simulation in the selected experiment.')
    parser_stdout.add_argument('--succeeded', action='store_true', help='Get the stdout for the first succeeded simulation in the selected experiment.')
    parser_stdout.set_defaults(func=func)


# 'dtk analyze' options
def populate_analyze_arguments(subparsers, func):
    parser_analyze = subparsers.add_parser('analyze', help='Analyze finished simulations in experiment according to analyzers.')
    parser_analyze.add_argument('-bn', '--batch_name', dest='batch_name', default=None, help='Use Batch Name for analyze.')
    parser_analyze.add_argument('-id', '--id', dest='itemids', default=None, nargs='+', help='IDs of the items to analyze (can be suites, batches, experiments, simulations)')
    parser_analyze.add_argument('-a', '--analyzer', default=None, help='Python script or builtin analyzer name for custom analysis of simulations (see dtk analyze-list).')
    parser_analyze.add_argument('-f', '--force', action='store_true', help='Force analyzer to run even if jobs are not all finished.')
    parser_analyze.set_defaults(func=func)


# 'dtk create_batch' options
def populate_createbatch_arguments(subparsers, func):
    parser_createbatch = subparsers.add_parser('create_batch', help='Create a Batch for later use in Analyze.')
    parser_createbatch.add_argument('-bn', '--batch_name', dest='batch_name', default=None, help='Use Batch Name.')
    parser_createbatch.add_argument('-id', '--id', dest='itemids', default=None, nargs='+', help='IDs of the items to analyze (can be suites, batches, experiments)')
    parser_createbatch.set_defaults(func=func)


# 'dtk list_batch' options
def populate_listbatch_arguments(subparsers, func):
    parser_listbatch=subparsers.add_parser('list_batch', help='Report recent 20 list of batches in Batch.')
    parser_listbatch.add_argument('-id', dest='id_or_name', help='Batch id or name.')
    parser_listbatch.set_defaults(func=func)


# 'dtk clean_batch' options
def populate_cleanbatch_arguments(subparsers, func):
    parser_cleanbatch = subparsers.add_parser('clean_batch', help='Remove all empty Batches.')
    parser_cleanbatch.set_defaults(func=func)


# 'dtk clear_batch' options
def populate_clearbatch_arguments(subparsers, func):
    parser_clearbatch = subparsers.add_parser('clear_batch', help='Remove all associated experiments from Batch given id or name.')
    parser_clearbatch.add_argument('-id', dest='id_or_name', required=True, help='Batch id or name.')
    parser_clearbatch.set_defaults(func=func)


# 'dtk delete_batch' options
def populate_deletebatch_arguments(subparsers, func):
    parser_deletebatch = subparsers.add_parser('delete_batch', help='Delete all Batches or Batch with given Batch ID/name.')
    parser_deletebatch.add_argument('-id', '--batch_id', dest='batch_id', help='Batch ID.')
    parser_deletebatch.set_defaults(func=func)


# 'dtk analyze-list` options
def populate_analyzer_list_arguments(subparsers, func):
    parser_analyze_list = subparsers.add_parser('analyze-list', help='List the available builtin analyzers.')
    parser_analyze_list.set_defaults(func=func)


# 'dtk sync' options
def populate_sync_arguments(subparsers, func):
    parser_sync = subparsers.add_parser('sync', help='Synchronize the COMPS database with the local database.')
    parser_sync.add_argument('-d', '--days',  help='Limit the sync to a certain number of days back', dest='days')
    parser_sync.add_argument('-id', help='Sync a specific experiment from COMPS.', dest='exp_id')
    parser_sync.add_argument('-n', '--name', help='Sync a specific experiment from COMPS (use %% for wildcard character).', dest='exp_name')
    parser_sync.add_argument('-u', '--user', help='Sync experiments belonging to a particular user', dest='user')
    parser_sync.set_defaults(func=func)


# 'dtk version' options
def populate_version_arguments(subparsers, func):
    parser_version = subparsers.add_parser('version', help='Display the current dtk-tools version.')
    parser_version.set_defaults(func=func)


# 'dtk test' options
def populate_test_arguments(subparsers, func):
    parser_test = subparsers.add_parser('test', help='Launch the nosetests on the test folder.')
    parser_test.set_defaults(func=func)


# 'dtk log' options
def populate_log_arguments(subparsers, func):
    parser_log = subparsers.add_parser('log', help="Allow to query and export the logs.")
    parser_log.add_argument('-l', '--level', help="Only display logs for a certain level and above (DEBUG,INFO,ERROR)", dest="level", default="DEBUG")
    parser_log.add_argument('-m', '--module', help="Only display logs for a given module.", dest="module", nargs='+')
    parser_log.add_argument('-n', '--number', help="Limit the number of entries returned (default is 100).", dest="number", default=100)
    parser_log.add_argument('-e', '--export', help="Export the log to the given file.", dest="export")
    parser_log.add_argument('-c', '--complete', help="Export the complete log to a CSV file (dtk_tools_log.csv).", action='store_true')
    parser_log.set_defaults(func=func)


# 'dtk list_packages' options
def populate_list_packages_arguments(subparsers, func):
    parser_list_packages = subparsers.add_parser('list_packages', help="List the packages available to get_package command.")
    parser_list_packages.set_defaults(func=func)


# 'dtk list_package_versions' options
def populate_list_package_versions_arguments(subparsers, func):
    parser = subparsers.add_parser('list_package_versions', help="List the versions available for this particular package.")
    parser.add_argument('package_name', help='The package to list versions of.')
    parser.set_defaults(func=func)


# 'dtk get_package' options
def populate_get_package_arguments(subparsers, func):
    parser = subparsers.add_parser('get_package', help="Obtain the specified disease configuration package for use.")
    parser.add_argument('package_name', help='The package name to obtain.')
    parser.add_argument('-v', '--version', help="Obtain a specific package version (Default: latest)", dest="package_version",
                        default='latest')
    parser.set_defaults(func=func)


# 'dtk generate_classes' options
def populate_generate_classes_arguments(subparsers, func):
    parser_sync = subparsers.add_parser('generate_classes', help='Generate classes and enums from Eradication SChema.')
    parser_sync.add_argument('-e', help='Path to the Eradication EXE.', required=True, dest='exe')
    parser_sync.add_argument('-o', help='Path of the output folder', required=True, dest='output')
    parser_sync.add_argument('-d', '--debug', help="Obtain a specific package version (Default: latest)", dest="debug",
                             action="store_true")
    parser_sync.set_defaults(func=func)