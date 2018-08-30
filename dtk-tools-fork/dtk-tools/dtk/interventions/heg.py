from dtk.utils.Campaign.CampaignClass import *


def heg_release(cb, released_number, num_repetitions = 13, timesteps_between_reps = 14):
    """
    Homing endonucleouse genetic intervention.

    :param cb: The :py:class:`DTKConfigBuilder <dtk.utils.core.DTKConfigBuilder>` object
    :return: Nothing
    """

    heg_release_event = CampaignEvent(
        Event_Coordinator_Config=StandardInterventionDistributionEventCoordinator(
            Intervention_Config=MosquitoRelease(
                Released_Gender=MosquitoRelease_Released_Gender_Enum.VECTOR_MALE,
                Released_Sterility=MosquitoRelease_Released_Sterility_Enum.VECTOR_FERTILE,
                Released_Genetics={
                    "Pesticide_Resistance": "WILD",
                    "HEG": "FULL"
                },
                Mated_Genetics={
                    "Pesticide_Resistance": "NotMated",
                    "HEG": "NotMated"
                },
                Released_HEGs=["FULL", "NotMated"],
                Released_Number=released_number,
                Released_Species="arabiensis",
                Cost_To_Consumer=200,
                Cost_To_Consumer_Citation="Alphey et al Vector Borne Zoonotic Dis 2010 10 295 by projecting 1979 An albimanus cost in El Salvador to 2008 dollars"
            ),
            Number_Repetitions=num_repetitions,
            Timesteps_Between_Repetitions=timesteps_between_reps
        ),
        Event_Name="MosquitoRelease",
        Nodeset_Config=NodeSetNodeList(Node_List=[340461476]),
        Start_Day=365
    )

    cb.add_event(heg_release_event)
    
    return {'num_released': released_number,
            'num_repetitions': num_repetitions,
            'timesteps_btw_reps': timesteps_between_reps}