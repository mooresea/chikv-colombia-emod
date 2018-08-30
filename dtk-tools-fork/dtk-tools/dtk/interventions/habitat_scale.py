import copy
from dtk.utils.Campaign.utils.RawCampaignObject import RawCampaignObject


def scale_larval_habitats(cb, df, start_day=0, repetitions=1, tsteps_btwn_repetitions=365,
                          node_property_restrictions=[]):
    """
    Reduce available larval habitat in a node-specific way.

    :param cb: The :py:class:`DTKConfigBuilder <dtk.utils.core.DTKConfigBuilder>` object
    :param df: dataframe containing habitat scale factors.
    Example 1: scale TEMPORARY_RAINFALL by 3-fold for all nodes, all species.
    df = pd.DataFrame( { 'TEMPORARY_RAINFALL': [3],
                     })
    Example 2: scale TEMPORARY_RAINFALL by 3-fold for all nodes, arabiensis only.
    df = pd.DataFrame( { 'TEMPORARY_RAINFALL.arabiensis': [3],
                     })
    Example 3: scale differently by node ID.
    df = pd.DataFrame( { 'NodeID' : [0, 1, 2, 3, 4],
                         'CONSTANT': [1, 0, 1, 1, 1],
                         'TEMPORARY_RAINFALL': [1, 1, 0, 1, 0],
                         })
    Example 4: scale differently by both node ID and species.
    df = pd.DataFrame( { 'NodeID' : [0, 1, 2, 3, 4],
                     'CONSTANT.arabiensis': [1, 0, 1, 1, 1],
                     'TEMPORARY_RAINFALL.arabiensis': [1, 1, 0, 1, 0],
                     'CONSTANT.funestus': [1, 0, 1, 1, 1]
                     })
    Example 5: scale some habitats by species and others same for all species.
    df = pd.DataFrame( { 'NodeID' : [0, 1, 2, 3, 4],
                         'CONSTANT.arabiensis': [1, 0, 1, 1, 1],
                         'TEMPORARY_RAINFALL.arabiensis': [1, 1, 0, 1, 0],
                         'CONSTANT.funestus': [1, 0, 1, 1, 1],
                         'LINEAR_SPLINE': [1, 1, 0, 1, 0]
                         })
    Example 6: scale nodes at different dates.
    df = pd.DataFrame( { 'NodeID' : [0, 1, 2, 3, 4],
                         'CONSTANT': [1, 0, 1, 1, 1],
                         'TEMPORARY_RAINFALL': [1, 1, 0, 1, 0],
                         'Start_Day': [0, 30, 60, 65, 65],
                         })
    :param start_day: date of scaling habitats for all scaling actions specified in df.
    Only applied if no Start_Day column in df.
    :return: Nothing
    """

    if 'Start_Day' not in df.columns.values:
        df['Start_Day'] = start_day

    standard_columns = ['NodeID', 'Start_Day']
    habitat_cols = [x for x in df.columns.values if x not in standard_columns]
    habitat_names = list(set([x.split('.')[0] for x in habitat_cols]))
    by_species = any(['.' in x for x in df.columns.values if x not in standard_columns])
    by_node = True if 'NodeID' in df.columns.values else False

    for start_day, df_by_date in df.groupby('Start_Day'):

        for gn, gdf in df_by_date.groupby(habitat_cols):
            if not by_species:
                if len(habitat_names) == 1:
                    hab_scales = {habitat_cols[0]: gn}
                else:
                    hab_scales = {x: y for x, y in zip(habitat_cols, gn)}
            else:
                if len(habitat_names) == 1:
                    if len(habitat_cols) == 1:
                        hab, sp = habitat_cols[0].split('.')
                        hab_scales = {hab: {sp: gn}}
                    else:
                        hab = habitat_cols[0].split('.')[0]
                        species = [x.split('.')[1] for x in habitat_cols]
                        hab_scales = {hab: {sp: x for (sp, x) in zip(species, gn)}}
                else:
                    hab_scales = {}
                    for ih, hab in enumerate(habitat_names) :
                        if hab in habitat_cols:
                            hab_scales[hab] = gn[ih]
                        else:
                            h = [x for x in habitat_cols if hab in x]
                            vals = [gn[x] for x in range(len(habitat_cols)) if hab in habitat_cols[x]]
                            hab_scales[hab] = {x.split('.')[1]: y for x, y in zip(h, vals)}

            if by_node:
                node_cfg = {
                    "class": "NodeSetNodeList",
                    "Node_List": [int(x) for x in gdf['NodeID']]
                }
            else:
                node_cfg = {"class": "NodeSetAll"}

            add_habitat_reduction_event(cb, start_day=start_day, node_cfg=node_cfg, hab_scales=hab_scales,
                                        repetitions=repetitions, tsteps_btwn_repetitions=tsteps_btwn_repetitions,
                                        node_property_restrictions=node_property_restrictions)


def add_habitat_reduction_event(cb, start_day, node_cfg, hab_scales, repetitions, tsteps_btwn_repetitions,
                                node_property_restrictions):

    # A permanent node-specific scaling of larval habitat by habitat type
    habitat_reduction_event = {
        "class": "CampaignEvent",
        "Start_Day": start_day,
        "Event_Coordinator_Config": {
            "class": "StandardInterventionDistributionEventCoordinator",
            "Number_Repetitions": repetitions,
            "Timesteps_Between_Repetitions": tsteps_btwn_repetitions,
            "Node_Property_Restrictions": node_property_restrictions,
            "Intervention_Config": {
                "class": "ScaleLarvalHabitat"
            }
        }
    }

    habitat_reduction_event['Nodeset_Config'] = node_cfg
    habitat_reduction_event['Event_Coordinator_Config']['Intervention_Config']['Larval_Habitat_Multiplier'] = hab_scales

    cb.add_event(RawCampaignObject(habitat_reduction_event))
