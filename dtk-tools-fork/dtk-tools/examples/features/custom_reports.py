from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.study_sites import configure_site
from malaria.reports.MalariaReport import add_summary_report, add_immunity_report, add_survey_report

from simtools.SetupParser import SetupParser
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory

exp_name = 'CustomReports'
cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')
configure_site(cb, 'Namawala')

nyears = 2
cb.update_params({'Num_Cores': 1,
                  'Base_Population_Scale_Factor': 0.1,
                  'x_Temporary_Larval_Habitat': 0.05,
                  'Simulation_Duration': 365 * nyears})

add_summary_report(cb, description='Monthly', interval=30)
add_summary_report(cb, description='Annual', interval=365)
add_immunity_report(cb, start=365 * (nyears - 1), interval=365, nreports=1, description="FinalYearAverage")
add_survey_report(cb, survey_days=[100, 200], reporting_interval=21, nreports=1)

run_sim_args = {'config_builder': cb,
                'exp_name': exp_name}

# If you prefer running with `python custom_reports.py`, you will need the following block
if __name__ == "__main__":
    SetupParser.init()
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    # Wait for the simulations to be done
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())
