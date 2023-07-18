"""Auxiliary functions for mapping to check content and set defaults of CUDS calculation objects."""
import os
import yaml
from osp.core.cuds import Cuds
from osp.tools.io_functions import raise_warning
from osp.core.namespaces import emmo


def AMS_default_setting(root_cuds_object: Cuds, accuracy_level: str,
                        calculation_type: str, setting: str) -> str:
    """Set default setting of a Wrapper object according to a yaml dictionary.

    Add the default setting as CUDS object in the wrapper object.

    :param root_cuds_object: = CUDS wrapper to be checked.

    :param accuracy_level: Level of accuracy.

    :param calculation_type: Calculation type.

    :param setting: Setting to be set by defaults.

    :return srt: Label of the default setting. It also adds the setting
                 as CUDS object in the wrapper.
    """

    filedir = os.path.dirname(os.path.abspath(__file__))
    dictionary = os.path.join(
        filedir, '../dictionaries/defaults/defaults_AMS.yaml')
    with open(dictionary, 'r') as f:
        AMS_list = yaml.load(f, Loader=yaml.FullLoader)
        default_setting = (AMS_list['AMS'][str(accuracy_level)][str(
            calculation_type)][setting])

    # Set default ExchangeCorrelationFunctional:
    if setting == 'ExchangeCorrelationFunctional':
        list_of_XCs = {
            'B3LYP': emmo.B3LYP,
            'LDA': emmo.LDA,
            'VWN': emmo.VWN,
            'PBE': emmo.PBE
        }
        xc_functional = list_of_XCs[default_setting]()
        root_cuds_object.add(xc_functional, rel=emmo.hasInput)

        raise_warning(
            file=os.path.basename(__file__),
            function=AMS_default_setting.__name__,
            message="Defining ExchangeCorrelationFunctional default, from defaults_AMS.yaml",
            message2=default_setting)

    # Set default BasisSet:
    if setting == 'BasisSet':
        list_of_basis = {
            'SZ': emmo.SZ,
            'DZ': emmo.DZ,
            'DZP': emmo.DZP,
            'TZP': emmo.TZP,
            'TZ2P': emmo.TZ2P
        }
        basis_set = list_of_basis[default_setting]()
        root_cuds_object.add(basis_set, rel=emmo.hasInput)

        raise_warning(
            file=os.path.basename(__file__),
            function=AMS_default_setting.__name__,
            message="Defining BasisSet default, from defaults_AMS.yaml",
            message2=default_setting)

    return default_setting
