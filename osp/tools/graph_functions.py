"""
Functions to graph dependencies between CUDS objects."""

import os
from osp.core.cuds import Cuds
from osp.core.ontology.relationship import OntologyRelationship
from osp.core.utils import simple_search as search
from osp.tools.io_functions import raise_error
from osp.core.namespaces import emmo
from osp.core.utils.general import get_relationships_between
from typing import List
# from osp.core.utils import pretty_print


def graph_wrapper_dependencies(root_cuds_object: Cuds) -> dict:
    """
    Check the Simulation/Calculation dependencies in a Wrapper object.

        :param root_cuds_object: Wrapper CUDS object.

        :return: dict containing the CUDS Calculations listed by relationships.
    """

    # Find Simulation objects in cuba.Wrapper:
    search_simulation = \
        search.find_cuds_object(criterion=lambda x:
                                x.is_a(oclass=emmo.Simulation),
                                root=root_cuds_object,
                                rel=emmo.hasPart,
                                find_all=True,
                                max_depth=1)

    # Zero Simulations found:
    if len(search_simulation) == 0:

        wrapper_dependencies = {"Calculation": [Cuds],
                                "Relationship": [OntologyRelationship]}

        [calculation_list, relationship_list] = graph_calculation_dependencies(
                                   root_cuds_object=root_cuds_object,
                                   simulation=False)

        wrapper_dependencies["Calculation"] = calculation_list
        wrapper_dependencies["Relationship"] = relationship_list

    # One Simulation found:
    elif len(search_simulation) == 1:
        # Instantiate the dictionary that will contain the dependencies:

        wrapper_dependencies = {"Simulation": {
                                "Calculations": [Cuds],
                                "Relationships": [OntologyRelationship]}}

        [calculation_list, relationship_list] = graph_calculation_dependencies(
                            root_cuds_object=search_simulation[0],
                            simulation=True)
        wrapper_dependencies["Simulation"]["Calculations"] = \
            calculation_list
        wrapper_dependencies["Simulation"]["Relationships"] = \
            relationship_list

    # More than one Simulation found:
    else:
        raise_error(file=os.path.basename(__file__),
                    function=graph_wrapper_dependencies.__name__,
                    type='NameError',
                    message="More than one Simulation have been found"
                            " in the Wrapper object.")
    return wrapper_dependencies


def graph_calculation_dependencies(root_cuds_object: Cuds,
                                   simulation: bool = False) -> tuple:
    """
    Graph dependencies of Calculations in a CUDS object.

        :param root_cuds_object: Wrapper CUDS object.

        :param simulation: boolean, is there and emmo.Simulation object?

        :return: tuple of list containing the Calculations CUDS
                and their relationships with the rood_cuds_object.
    """

    calculation_list = []
    relationship_list = []
    calculations_types = emmo.Calculation.subclasses
    if simulation:
        first_calculation = \
            search.find_relationships(find_rel=emmo.INVERSE_OF_hasTemporalFirst,
                                      root=root_cuds_object,
                                      consider_rel=emmo.EMMORelation)
    else:
        first_calculation = []

    if (first_calculation):
        calculation_list.append(first_calculation[0])
        relationship_list.append(emmo.hasTemporalFirst)

        next_calculation = \
            search.find_relationships(find_rel=emmo.INVERSE_OF_hasTemporalNext,
                                      root=root_cuds_object,
                                      consider_rel=emmo.hasTemporalNext)

        if next_calculation:
            calculation_list.append(next_calculation[0])
            relationship_list.append(emmo.hasTemporalNext)
            last_calculation = \
                search.find_relationships(find_rel=emmo.INVERSE_OF_hasTemporalLast,
                                          root=root_cuds_object,
                                          consider_rel=emmo.EMMORelation)
            if last_calculation:
                calculation_list.append(last_calculation[0])
                relationship_list.append(emmo.hasTemporalLast)
            else:
                raise_error(file=os.path.basename(__file__),
                            function=graph_calculation_dependencies.__name__,
                            type='NameError',
                            message="emmo.hasTemporalNext calculation without"
                                    " an emmo.hasTemporalLast.")
        else:
            raise_error(file=os.path.basename(__file__),
                        function=graph_calculation_dependencies.__name__,
                        type='NameError',
                        message="emmo.hasTemporalFirst calculation without"
                                " an emmo.hasTemporalNext.")
    else:
        search_calculation = \
            search.find_cuds_object(criterion=lambda x:
                                    x.oclass in calculations_types,
                                    root=root_cuds_object,
                                    rel=emmo.hasPart,
                                    find_all=True,
                                    max_depth=1)
        relationship_list.append(
            get_relationships_between(root_cuds_object, search_calculation[0]))

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

    return calculation_list, relationship_list


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
