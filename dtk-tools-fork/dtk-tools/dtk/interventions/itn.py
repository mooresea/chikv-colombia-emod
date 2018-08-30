import copy
from dtk.interventions.triggered_campaign_delay_event import triggered_campaign_delay_event
from dtk.utils.Campaign.CampaignClass import *


itn_bednet = SimpleBednet(Bednet_Type='ITN',
                          Killing_Config=WaningEffectExponential(Initial_Effect=0.6, Decay_Time_Constant=1460),
                          Blocking_Config=WaningEffectExponential(Initial_Effect=0.9, Decay_Time_Constant=730),
                          Usage_Config=WaningEffectRandomBox(Expected_Discard_Time=3650, Initial_Effect=1.0),
                          Cost_To_Consumer=3.75
                          )

receiving_itn_event = BroadcastEvent(Broadcast_Event='Received_ITN')


def add_ITN(config_builder, start, coverage_by_ages, waning={}, cost=None, nodeIDs=[], node_property_restrictions=[],
            ind_property_restrictions=[], triggered_campaign_delay=0, trigger_condition_list=[], listening_duration=1):
    """
     Add an ITN intervention to the config_builder passed.
    :param config_builder: The :py:class:`DTKConfigBuilder <dtk.utils.core.DTKConfigBuilder>` holding the campaign that will receive the ITN event
    :param start: The start day of the bednet distribution
    :param coverage_by_ages: a list of dictionaries defining the coverage per age group
        [{"coverage":1,"min": 1, "max": 10},{"coverage":1,"min": 11, "max": 50},{ "coverage":0.5, "birth":"birth", "duration":34}]
        birthtriggered(in coverage_by_age) and triggered_condition_list are mututally exclusive. "birth" option will be ingnored if you're
        using trigger_condition_list
    :param waning: a dictionary defining the durability of the nets. if empty the default decay profile will be used.
    For example, update usage duration to 180 days as waning={'Usage_Config' : {"Expected_Discard_Time": 180}}
    :param cost: Set the ``Cost_To_Consumer`` parameter
    :param nodeIDs: If empty, all nodes will get the intervention. If set, only the nodeIDs specified will receive the intervention.
    :param node_property_restrictions: restricts itn based on list of node properties in format [{"Place":"RURAL"}, {"ByALake":"Yes, "LovelyWeather":"No}]
    :param ind_property_restrictions: Restricts itn based on list of individual properties in format [{"BitingRisk":"High", "IsCool":"Yes}, {"IsRich": "Yes"}]
    :param triggered_campaign_delay: how many time steps after receiving the trigger will the campaign start.
    Eligibility of people or nodes for campaign is evaluated on the start day, not the triggered day.
    :param trigger_condition_list: when not empty,  the start day is the day to start listening for the trigger conditions listed, distributing the spraying
        when the trigger is heard. This does not distribute the BirthTriggered intervention.
    :param listening_duration: how long the distributed event will listen for the trigger for, default is 1, which is indefinitely
    :return: Nothing
    """

    if waning:
        for w, w_config in waning.items():
            setattr(itn_bednet, w, w_config)

    itn_bednet.Cost_To_Consumer = cost

    itn_bednet_w_event = MultiInterventionDistributor(Intervention_List=[itn_bednet, receiving_itn_event])

    # Assign node IDs #
    # Defaults to all nodes unless a node set is specified
    if not nodeIDs:
        nodeset_config = NodeSetAll()
    else:
        nodeset_config = NodeSetNodeList(Node_List=nodeIDs)

    if triggered_campaign_delay:
        trigger_condition_list = [str(triggered_campaign_delay_event(config_builder,
                                                                     start,  nodeIDs,
                                                                     triggered_campaign_delay,
                                                                     trigger_condition_list,
                                                                     listening_duration))]

    for coverage_by_age in coverage_by_ages:
        if trigger_condition_list:
            if not 'birth' in coverage_by_age.keys():
                intervention_config = NodeLevelHealthTriggeredIV(
                    Trigger_Condition_List=trigger_condition_list,
                    Duration=listening_duration,
                    Demographic_Coverage=coverage_by_age["coverage"],
                    Target_Residents_Only=True,
                    Actual_IndividualIntervention_Config=itn_bednet_w_event  # itn_bednet
                )

                ITN_event = CampaignEvent(Start_Day=int(start),
                                          Nodeset_Config=nodeset_config,
                                          Event_Coordinator_Config=StandardInterventionDistributionEventCoordinator(
                                              Intervention_Config=intervention_config)
                                          )

                if all([k in coverage_by_age.keys() for k in ['min', 'max']]):
                    ITN_event_e_i = ITN_event.Event_Coordinator_Config.Intervention_Config
                    ITN_event_e_i.Target_Demographic = StandardInterventionDistributionEventCoordinator_Target_Demographic_Enum.ExplicitAgeRanges
                    ITN_event_e_i.Target_Age_Min = coverage_by_age["min"]
                    ITN_event_e_i.Target_Age_Max = coverage_by_age["max"]

                if ind_property_restrictions:
                    ITN_event_e_i.Property_Restrictions_Within_Node = ind_property_restrictions

                if node_property_restrictions:
                    ITN_event_e_i.Node_Property_Restrictions = node_property_restrictions

        else:
            event_coordinator_config = StandardInterventionDistributionEventCoordinator(
                Node_Property_Restrictions=[],
                Target_Residents_Only=1,
                Demographic_Coverage=coverage_by_age["coverage"],
                Intervention_Config=itn_bednet_w_event              # itn_bednet
            )
            ITN_event = CampaignEvent(Start_Day=int(start),
                                      Nodeset_Config=nodeset_config,
                                      Event_Coordinator_Config=event_coordinator_config
                                      )

            if node_property_restrictions:
                ITN_event.Event_Coordinator_Config.Node_Property_Restrictions.extend(node_property_restrictions)

            if all([k in coverage_by_age.keys() for k in ['min', 'max']]):
                ITN_event.Event_Coordinator_Config.Target_Demographic = StandardInterventionDistributionEventCoordinator_Target_Demographic_Enum.ExplicitAgeRanges
                ITN_event.Event_Coordinator_Config.Target_Age_Min = coverage_by_age["min"]
                ITN_event.Event_Coordinator_Config.Target_Age_Max = coverage_by_age["max"]

            if 'birth' in coverage_by_age.keys() and coverage_by_age['birth']:
                birth_triggered_intervention = BirthTriggeredIV(
                    Duration=coverage_by_age.get('duration', -1),               # default to forever if  duration not specified
                    Demographic_Coverage=coverage_by_age["coverage"],
                    Actual_IndividualIntervention_Config=itn_bednet_w_event     # itn_bednet
                )

                ITN_event.Event_Coordinator_Config.Intervention_Config = birth_triggered_intervention
                del ITN_event.Event_Coordinator_Config.Demographic_Coverage
                del ITN_event.Event_Coordinator_Config.Target_Residents_Only

                if ind_property_restrictions:
                    ITN_event.Event_Coordinator_Config.Intervention_Config.Property_Restrictions_Within_Node = ind_property_restrictions

            elif ind_property_restrictions:
                ITN_event.Event_Coordinator_Config.Property_Restrictions_Within_Node = ind_property_restrictions

        config_builder.add_event(ITN_event)