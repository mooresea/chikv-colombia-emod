import os
import re


class CMSParser:
    """
    With current version it only supports the following key CMS components:
     - start-model
     - reaction
     - param
     - func
     - bool
     - observe
     - species
     - time-event
     - state-event
    """

    re_map = {
        'start-model': """\(\s*(?P<key>start-model)\s+(?P<quote>['"])(?P<name>\S+)(?P=quote)\s*\)""",
        'param': r'\(\s*(?P<key>param)\s+(?P<name>\S+)\s+(?P<value>.*)\s*\)',
        'func': r'\(\s*(?P<key>func)\s+(?P<name>\S+)\s+(?P<func>.*)\s*\)',
        'bool': r'\(\s*(?P<key>bool)\s+(?P<name>\S+)\s+(?P<expr>.*)\s*\)',
        'observe': r'\(\s*(?P<key>observe)\s+(?P<label>\S+)\s+(?P<func>.*)\s*\)',
        'species1': r'\(\s*(?P<key>species)\s+(?P<name>\S+)\s*\)',
        'species2': r'\(\s*(?P<key>species)\s+(?P<name>\S+)\s+(?P<value>.+)\s*\)',
        'reaction': r'\(\s*(?P<key>reaction)\s+(?P<name>\S+)\s+(?P<input>\(.*?\))\s+(?P<output>\(.*?\))\s+(?P<func>.*)\s*\)',
        'time-event1': r'\(\s*(?P<key>time-event)\s+(?P<name>\S+)\s+(?P<time>\S+)\s+\((?P<pairs>\(.*\))\s*\)',
        'time-event2': r'\(\s*(?P<key>time-event)\s+(?P<name>\S+)\s+(?P<time>\S+)\s+(?P<iteration>\S+)\s+\((?P<pairs>\(.*\))\s*\)',
        'state-event': r'\(\s*(?P<key>state-event)\s+(?P<name>\S+)\s+(?P<predicate>.*?)\s+\((?P<pairs>.*)\)\s*\)'
    }

    reg_pair = r'(\s*\(\s*(.*)\s*\))+'

    @staticmethod
    def parse_model_from_file(model_file, cb):
        """
        parse model and build up cb
        """
        # normalize the paths
        model_file = os.path.abspath(os.path.normpath(model_file))

        # read the model file
        model = open(model_file, 'r').read()

        # clean model file
        model = CMSParser.clean_model(model)

        # parse the model
        CMSParser.parse_model(model, cb)

    @staticmethod
    def parse_model(model, cb):
        """
        parse model and populate cb properties
        """
        for exp in CMSParser.re_map.values():
            reg = re.compile(exp)

            matches = reg.finditer(model)
            for match in matches:
                term = match.groups()
                if term[0] == 'start-model':
                    cb.start_model_name = term[2]
                elif term[0] == 'species':
                    if len(term) == 3:
                        cb.add_species(term[1], term[2])
                    else:
                        cb.add_species(term[1])
                elif term[0] == 'param':
                    cb.add_param(term[1], term[2])
                elif term[0] == 'func':
                    cb.add_func(term[1], term[2])
                elif term[0] == 'bool':
                    cb.add_bool(term[1], term[2])
                elif term[0] == 'observe':
                    cb.add_observe(term[1], term[2])
                elif term[0] == 'reaction':
                    cb.add_reaction(term[1], term[2], term[3], term[4])
                elif term[0] == 'time-event':
                    if len(term) == 5:
                        cb.add_time_event(term[1], term[2], term[3], *CMSParser.parse_pairs(term[4]))
                    else:
                        cb.add_time_event(term[1], term[2], None, *CMSParser.parse_pairs(term[3]))
                elif term[0] == 'state-event':
                    cb.add_state_event(term[1], term[2], *CMSParser.parse_pairs(term[3]))
                else:
                    raise Exception('{} is not supported at moment.'.format([0]))

            # clear parsed text (workaround to distinguish two species cases)
            model = CMSParser.clean_model(model, exp)

    @staticmethod
    def parse_pairs(line):
        """
        parse pair string and return as a list
        """
        reg = re.compile(CMSParser.reg_pair)
        matches = reg.finditer(line)

        pairs = []
        for match in matches:
            pairs.extend(list(match.groups()))

        return pairs

    @staticmethod
    def clean_model(model, expr=r';.*'):
        """
        default: remove comments
        """
        return re.sub(expr, '', model)
