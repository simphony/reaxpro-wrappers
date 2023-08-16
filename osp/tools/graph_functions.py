"""
Functions to graph dependencies between CUDS objects."""

import os
from osp.core.ontology import OntologyClass
from osp.core.cuds import Cuds
from osp.core.utils import simple_search as search
from osp.tools.io_functions import raise_error
from osp.core.namespaces import emmo, cuba
from typing import List
# from osp.core.utils import pretty_print


def graph_wrapper_dependencies(root_cuds_object: Cuds, oclass=emmo.AtomisticCalculation) -> dict:
    """
    Check the Simulation/Calculation dependencies in a Wrapper object.

        :param root_cuds_object: Wrapper CUDS object.

        :return: dict containing the CUDS Calculations listed by relationships.
    """

    # Find Simulation objects in cuba.Wrapper:
    search_simulation = \
        search.find_cuds_object(criterion=lambda x:
                                x.is_a(oclass=emmo.Simulation) or x.is_a(emmo.Workflow),
                                root=root_cuds_object,
                                rel=cuba.relationship,
                                find_all=True,
                                max_depth=1)

    # Zero Simulations found:
    if len(search_simulation) == 0:


        calculation_list = graph_calculation_dependencies(root_cuds_object, oclass)

        wrapper_dependencies = {"Calculation": calculation_list}

    # One Simulation found:
    elif len(search_simulation) == 1:
        # Instantiate the dictionary that will contain the dependencies:

        calculation_list = graph_calculation_dependencies(search_simulation[0], oclass)
        wrapper_dependencies = {
                                    "Simulation": {
                                        "Calculations": calculation_list
                                    }
                               }


    # More than one Simulation found:
    else:
        raise_error(file=os.path.basename(__file__),
                    function=graph_wrapper_dependencies.__name__,
                    type='NameError',
                    message="More than one Simulation have been found"
                            " in the Wrapper object.")
    return wrapper_dependencies


def graph_calculation_dependencies(root_cuds_object: Cuds, oclass: OntologyClass) -> list:
    """
    Graph dependencies of Calculations in a CUDS object.

        :param root_cuds_object: Wrapper CUDS object.

        :param simulation: boolean, is there and emmo.Simulation object?

        :return: list containing the Calculations CUDS
    """

    calculation_list = []
    if root_cuds_object.is_a(emmo.Simulation) or root_cuds_object.is_a(emmo.Workflow):
        first_calculation = \
            search.find_cuds_object(criterion=lambda x:
                                    x.is_a(oclass),
                                    root=root_cuds_object,
                                    rel=cuba.relationship,
                                    find_all=True,
                                    max_depth=1)
    else:
        first_calculation = []

    if first_calculation:
        current = first_calculation

        while current:
             current = current.pop()
             calculation_list.append(current)
             current = current.get(oclass=oclass, rel=emmo.hasSpatialNext) \
                or current.get(oclass=emmo.PostProcessing, rel=emmo.hasSpatialNext)

    else:
        search_calculation = \
            search.find_cuds_object(criterion=lambda x:
                                    x.is_a(oclass),
                                    root=root_cuds_object,
                                    rel=cuba.relationship,
                                    find_all=True,
                                    max_depth=1)
    
        if len(search_calculation) == 0:
            raise_error(file=os.path.basename(__file__),
                        function=graph_calculation_dependencies.__name__,
                        type='NameError',
                        message="No Calculation found in CUDS object.")

        elif len(search_calculation) == 1:
            calculation_list.append(search_calculation[0])

        elif len(search_calculation) >= 1:
            calculation_list = solve_calculation_dependencies(
                                    cuds_list=search_calculation)
    return calculation_list


def solve_calculation_dependencies(cuds_list: List[Cuds]) -> list:
    """
    Set order of execution in a cuds_list of Calculations according to engine.

        :param cuds_list: List of Calculations CUDS object.

        :return: list of Calculation uids in order of execution.
    """

    solved_list = []
    # TODO A very simple schema to solve dependencies, to be improved:

    for cl in cuds_list:
        if cl.oclass == emmo.WavefunctionOptimization:
            solved_list.append(cl.uid)

        if cl.oclass == emmo.GeometryOptimization:
            solved_list.append(cl.uid)

        if cl.oclass == emmo.StationaryPointCalculation:
            solved_list.append(cl.uid)

        if cl.oclass == emmo.VibrationalFrequencyCalculation:
            solved_list.append(cl.uid)

        if cl.oclass == emmo.PotentialEnergySurfaceScan:
            solved_list.append(cl.uid)

        if cl.oclass == emmo.Calculation:
            solved_list.append(cl.uid)

    return solved_list
