"""Functions to map data between PLAMS and the wrapper object."""
import os
import shutil
import tempfile
import scm.pyzacros as pz
from scm.plams import MultiJob, AMSJob
from scm.plams import Settings as PlamsSettings
from scm.plams import Molecule as PlamsMolecule
from scm.plams import Atom as PlamsAtom
from arcp import arcp_random
from osp.tools.io_functions import raise_error, raise_warning
from osp.tools.io_functions import read_mechanism, read_cluster_expansion, read_molecule
from osp.tools.set_functions import AMS_default_setting
from osp.models.utils.general import get_upload
from osp.core.utils import simple_search as search
from osp.core.ontology.oclass import OntologyClass
from osp.core.cuds import Cuds
from osp.core.namespaces import emmo, crystallography, cuba
from osp.core.utils import pretty_print
from typing import List
from uuid import UUID

def map_function(self, root_cuds_object: Cuds, engine=None) -> tuple:
    """
    Map a Calculations CUDS object to a engine-understandable input.

    :param root_cuds_object: Calculation CUDS object.

    :param engine: Engine provided when initializing the wrapper object, it can
                   be of different type.

    :param return: Returning tuple of the type depending on the engine chosen.
    """

    if isinstance(self.engine, MultiJob):

        software = emmo.AMS()
        root_cuds_object.add(software, rel=emmo.hasModellingSoftware)

        search_landscape = \
            search.find_cuds_objects_by_oclass(emmo.EnergyLandscape, root_cuds_object, emmo.hasPart)

        if search_landscape:
            plams_molecule = map_PLAMSLandscape(search_landscape[0])

        else:
            plams_molecule = map_PLAMSMolecule(root_cuds_object)

        path = os.path.join(self.workdir, self.jobname)
        plams_settings = map_PLAMSSettings(path, root_cuds_object)

        return (plams_molecule, plams_settings)

    elif self.engine == "Zacros":

        software = emmo.Zacros()
        root_cuds_object.add(software, rel=emmo.hasModellingSoftware)

        pz_settings = map_PyZacrosSettings(self.calculation)

        if hasattr(self, 'mechanism'):
            # if mechanism is read from CUDS graph
            pz_mechanism = map_PyZacrosMechanism(self.mechanism)
        else:
            # otherwise, use default location
            input_job = pz.ZacrosJob.load_external(self.input_path)
            pz_mechanism = input_job.mechanism

        if hasattr(self, 'lattice'):
            pz_lattice = pz.Lattice(fileName=self.lattice)
        else:
            pz_lattice = input_job.lattice

        if hasattr(self, 'cluster'):
            pz_cluster_expansions = map_PyZacrosClusterExpansion(self.cluster)
        else:
            pz_cluster_expansions = input_job.cluster_expansion

        return (pz_settings, pz_lattice, pz_mechanism, pz_cluster_expansions)

    else:
        # Raise error if the engine is other than AMS (to be extended of course!):
        raise_error(file=os.path.basename(__file__), function=map_function.__name__,
                    type='NameError', message='Engine not implemented yet.')
    return


def map_PLAMSMolecule(root_cuds_object: Cuds) -> PlamsMolecule:
    """
    Map ontology-defined molecule to PlamsMolecule() object.

    :param root_cuds_object: Calculation CUDS object.

    :return PlamsMolecule:
    """

    dict_molecule = map_molecule(root_cuds_object)
    dict_lattice = map_lattice(root_cuds_object)
    syntactic_molecule = PlamsMolecule()

    natoms = len(dict_molecule["atom_symbol"])

    if "region" in dict_molecule:

        for ij in range(natoms):
            syntactic_molecule.add_atom(
                PlamsAtom(
                    symbol=dict_molecule["atom_symbol"][ij],
                    coords=dict_molecule["atom_coordinates"][ij]))

        for index, ij in enumerate(syntactic_molecule.atoms):
            ij.properties.suffix = "region=" + dict_molecule["region"][index]

    else:

        for ij in range(natoms):
            syntactic_molecule.add_atom(PlamsAtom(
                                    symbol=dict_molecule["atom_symbol"][ij],
                                    coords=dict_molecule["atom_coordinates"][ij]))

    if "molecular_charge" in dict_molecule:
        syntactic_molecule.properties.charge = dict_molecule["molecular_charge"][0]

    if bool(dict_lattice):

        for ij in dict_lattice['vector_coordinates']:
            syntactic_molecule.lattice.append(ij)

    return syntactic_molecule


def map_PLAMSLandscape(root_cuds_object: Cuds) -> dict:
    """
    Map ontology-defined EnergyLandscape to a dict {PLAMSMolecule}.

    :param root_cuds_object: EnergyLandscape CUDS object.

    :return {PlamsMolecule}:
    """

    landscape_dict = {}
    search_reaction_path = find_reaction_path(root_cuds_object)

    if search_reaction_path:

        # Generate list of labels for single-molecule states:
        path_label_list = []

        for index, reaction_path in enumerate(search_reaction_path):

            path_label = "state" + str(1 + index)

            search_molecule = search.find_cuds_object(
                criterion=lambda x:
                x.is_a(oclass=emmo.MolecularGeometry), root=reaction_path, rel=emmo.hasPart,
                find_all=True, max_depth=1)
            if search_molecule:

                if len(search_molecule) > 1:
                    raise_error(
                        file=os.path.basename(__file__), function=map_PLAMSLandscape.__name__,
                        type='NameError', message='More than one Molecule '
                        'defined in ReactionPathway CUDS object.')

                path_tuple = (path_label, search_molecule[0].uid)
                path_label_list.append(path_tuple)
                landscape_dict[path_label] = map_PLAMSMolecule(search_molecule[0])

            else:

                search_reactants = \
                    search.find_cuds_objects_by_oclass(
                                       emmo.ChemicalReactionEquationReactant,
                                       reaction_path, emmo.hasPart)
                search_products = \
                    search.find_cuds_objects_by_oclass(
                                       emmo.ChemicalReactionEquationProduct,
                                       reaction_path, emmo.hasPart)
                search_TS = \
                    search.find_cuds_objects_by_oclass(
                                       emmo.TransitionStateGeometry,
                                       reaction_path, emmo.hasPart)

                if search_reactants:

                    for reactant in search_reactants:
                        search_molecule = \
                           search.find_cuds_objects_by_oclass(
                                                emmo.MolecularGeometry,
                                                reactant, emmo.hasPart)

                        if len(search_molecule) > 1:
                            raise_error(
                                file=os.path.basename(__file__),
                                function=map_PLAMSLandscape.__name__,
                                type='NameError',
                                message='More than one Molecule defined in '
                                'ChemicalReactionEquationReactant CUDS object.')

                        for label in path_label_list:
                            if search_molecule[0].uid == label[1]:
                                reactant_label = label[0]

                if search_products:

                    for product in search_products:
                        search_molecule = \
                           search.find_cuds_objects_by_oclass(
                                                emmo.MolecularGeometry,
                                                product, emmo.hasPart)
                        if len(search_molecule) > 1:
                            raise_error(file=os.path.basename(__file__),
                                        function=map_PLAMSLandscape.__name__,
                                        type='NameError',
                                        message='More than one Molecule defined in '
                                        'ChemicalReactionEquationProduct CUDS object.')

                        for label in path_label_list:

                            if search_molecule[0].uid == label[1]:
                                product_label = label[0]

                if search_TS:

                    if len(search_TS) > 1:
                        raise_error(file=os.path.basename(__file__),
                                    function=map_PLAMSLandscape.__name__,
                                    type='NameError',
                                    message='More than one TransitionStateGeometry defined '
                                            'in ReactionPathway CUDS object.')

                    search_molecule = \
                        search.find_cuds_objects_by_oclass(emmo.MolecularGeometry,
                                                           search_TS[0], emmo.hasPart)

                    if len(search_molecule) > 1:
                        raise_error(file=os.path.basename(__file__),
                                    function=map_PLAMSLandscape.__name__, type='NameError',
                                    message='More than one Molecule defined'
                                            ' in TransitionStateGeometry CUDS object.')
                    path_label = path_label + ' ts=T ' + 'reactant=' + \
                        reactant_label + ' product=' + product_label
                    landscape_dict[path_label] = map_PLAMSMolecule(search_molecule[0])

    else:
        raise_error(file=os.path.basename(__file__), function=map_molecule.__name__,
                    type='NameError', message='ReactionPathway CUDS object missed in the wrapper.')

    return landscape_dict


def map_molecule(root_cuds_object: Cuds) -> dict:
    """
    Map Molecule from a Calculation CUDS object.

    It returns a dictionary containing the geometry information:
    {
     "atom_symbol": [atom_label_1, atom_label_2, ..., atom_label_n]
     "atom_coordinates":
                    [[coordinate_x_1, coordinate_y_1, coordinate_z_1, (region)],
                    [coordinate_x_2, coordinate_y_2, coordinate_z_2, (region)],
                     ...
                    [coordinate_x_2, coordinate_y_n, coordinate_z_n, (region)]]
    }

    :param root_cuds_object: CUDS root_cuds_object to be checked.

    :return str: If a Molecule object is present, return dictionary with
                 atom labels and positions. Otherwise, raise an error.
    """

    search_molecule = \
        search.find_cuds_objects_by_oclass(emmo.MolecularGeometry, root_cuds_object, emmo.hasPart)

    if len(search_molecule) > 1:
        raise_error(file=os.path.basename(__file__),
                    function=map_function.__name__,
                    type='ValueError',
                    message='More than one emmo.MolecularGeometry defined in the Wrapper object.')
    elif not search_molecule:
        raise_error(file=os.path.basename(__file__), function=map_molecule.__name__,
                    type='ValueError', message='Molecule CUDS object missed in the wrapper.')

    molecule_data = {}
    # Looking for atoms:
    first = search_molecule.pop().get(rel=emmo.hasSpatialFirst)
    if len(first) > 1:
        raise_error(file=os.path.basename(__file__),
                    function=map_molecule.__name__,
                    type='ValueError',
                    message='More than one first spatial atom found in molecule.')
    elif len(first) == 0:
        raise_error(file=os.path.basename(__file__),
                    function=map_molecule.__name__,
                    type='ValueError',
                    message='First spatial atom in molecule not found.')
    else:
        atom_list = []
        current = first
        while current:
            current = current.pop()
            atom_list.append(current)
            current = current.get(rel=emmo.hasSpatialNext)
    # Looking for charges:
    search_charge =  \
        search.find_cuds_objects_by_oclass(emmo.ElectricCharge, root_cuds_object,
                                            emmo.hasProperty)

    # Is there any charge?
    if search_charge:
        molecular_charge = search_charge[0].get(oclass=emmo.Integer,
                                                rel=emmo.hasQuantityValue)[0].hasNumericalData

    # For each of the atoms, we add labels and positions:
    for count, iatom in enumerate(atom_list):

        search_label = \
            search.find_cuds_objects_by_oclass(emmo.ChemicalElement, iatom, emmo.hasPart)

        search_position_vec = \
            search.find_cuds_objects_by_oclass(emmo.PositionVector, iatom, emmo.hasPart)

        search_region = \
            search.find_cuds_objects_by_oclass(emmo.Region, iatom, emmo.hasPart)

        # Is there any region defined? For EON PES exploration.
        if search_region:
            region = search_region[0].hasSymbolData

        # Atom coordinates:
        coordinate_x = \
            search_position_vec[0].get(oclass=emmo.Real, rel=emmo.hasSpatialFirst)
        coordinate_y = \
            search_position_vec[0].get(oclass=emmo.Real, rel=emmo.hasSpatialNext)
        coordinate_z = \
            search_position_vec[0].get(oclass=emmo.Real, rel=emmo.hasSpatialLast)

        if count == 0:
            molecule_data = {
                'atom_symbol': [search_label[0].hasSymbolData],
                'atom_coordinates': [[
                    coordinate_x[0].hasNumericalData,
                    coordinate_y[0].hasNumericalData,
                    coordinate_z[0].hasNumericalData
                ]]
            }
            if search_region:
                molecule_data['region'] = [region]
            if search_charge:
                molecule_data['molecular_charge'] = [molecular_charge]
        else:
            molecule_data['atom_symbol'].\
                append(search_label[0].hasSymbolData)
            molecule_data['atom_coordinates'].\
                append([coordinate_x[0].hasNumericalData,
                        coordinate_y[0].hasNumericalData,
                        coordinate_z[0].hasNumericalData
                        ])
            if search_region:
                molecule_data['region'].append(region)

    return molecule_data


def map_lattice(root_cuds_object: Cuds) -> dict:
    """Map Lattice from a Calculation CUDS object.

    It returns a dictionary containing 3x3 lattice vector array
    {
     "lattice_vector":
                    [[coordinate_x_1, coordinate_y_1, coordinate_z_1],
                    [coordinate_x_2, coordinate_y_2, coordinate_z_2],
                    [coordinate_x_3, coordinate_y_3, coordinate_z_3]]
    }

    :param root_cuds_object: CUDS wrapper to be checked.

    :return str: If a Lattice object is present, return dictionary with
                 the lattice vectors. Otherwise, raise an error.
    """

    search_lattice = \
        search.find_cuds_objects_by_oclass(
                           crystallography.UnitCell,
                           root_cuds_object, emmo.hasPart)

    #  Raise error if two molecules, for some reason, are found:
    if len(search_lattice) > 1:
        raise_error(file=os.path.basename(__file__), function=map_molecule.__name__,
                    type='NameError', message='More than one crystallography.UnitCell'
                    ' defined in the Wrapper object.')

    if search_lattice:
        lattice_data = {}

        # Looking for atoms:
        first_vector = \
            search.find_relationships(find_rel=emmo.INVERSE_OF_hasSpatialFirst,
                                      root=search_lattice[0],
                                      consider_rel=emmo.hasSpatialFirst)
        next_vector = \
            search.find_relationships(find_rel=emmo.INVERSE_OF_hasSpatialNext,
                                      root=search_lattice[0],
                                      consider_rel=emmo.hasSpatialNext)

        last_vector = \
            search.find_relationships(find_rel=emmo.INVERSE_OF_hasSpatialLast,
                                      root=search_lattice[0],
                                      consider_rel=emmo.hasSpatialLast)

        vector_list = []
        if first_vector:
            vector_list.append(first_vector[0])
        if next_vector:
            vector_list.append(next_vector[0])
        if last_vector:
            vector_list.append(last_vector[0])

        # For each of them, we add labels and positions:
        for count, ivector in enumerate(vector_list):

            # Vector coordinates:
            coordinate_x = \
                ivector.get(oclass=emmo.Real, rel=emmo.hasSpatialFirst)
            coordinate_y = \
                ivector.get(oclass=emmo.Real, rel=emmo.hasSpatialNext)
            coordinate_z = \
                ivector.get(oclass=emmo.Real, rel=emmo.hasSpatialLast)
            if count == 0:
                lattice_data = {
                    'vector_coordinates': [(
                        float(coordinate_x[0].hasNumericalData),
                        float(coordinate_y[0].hasNumericalData),
                        float(coordinate_z[0].hasNumericalData))]
                                }
            else:
                lattice_data['vector_coordinates'].\
                    append((float(coordinate_x[0].hasNumericalData),
                            float(coordinate_y[0].hasNumericalData),
                            float(coordinate_z[0].hasNumericalData)))
    else:
        raise_warning(file=os.path.basename(__file__),
                      function=map_lattice.__name__,
                      message='UnitCell CUDS object missed in the wrapper.')
        return
    return lattice_data


def map_PLAMSSettings(workdir: str, root_cuds_object: Cuds) -> PlamsSettings:
    """
    Map ontology-defined calculation settings to PlamsSettings() object.

    :param root_cuds_object: Calculation CUDS object.

    :return PlamsSettings:
    """

    semantic_settings = {}
    syntactic_settings = PlamsSettings()

    # Map accuracy level:
    semantic_settings['AccuracyLevel'] = map_accuracy_level(root_cuds_object)

    # Map calculation type:
    semantic_settings['Calculation'] = map_calculation_type(root_cuds_object)

    # Map Model type:
    semantic_settings['Model'] = map_model_type(root_cuds_object)

    if semantic_settings['Calculation'] == "WavefunctionOptimization":
        syntactic_settings.input.AMS.task = "SinglePoint"

    elif semantic_settings['Calculation'] == "GeometryOptimization":
        syntactic_settings.input.AMS.task = "GeometryOptimization"

    elif semantic_settings['Calculation'] == "PotentialEnergySurfaceScan":
        syntactic_settings.input.AMS.task = "PESScan"

    elif semantic_settings['Calculation'] == "LandscapeRefinement":

        syntactic_settings.input.AMS.task = "PESExploration"
        syntactic_settings.input.AMS.GeometryOptimization.InitialHessian.Type = "Calculate"
        syntactic_settings.input.AMS.PESExploration.Job = "LandscapeRefinement"

        if semantic_settings['Model'] == "DFTB":
            syntactic_settings.input.DFTB = PlamsSettings()

    elif semantic_settings['Calculation'] == 'ProcessSearch':

        syntactic_settings.input.AMS.task = "PESExploration"

        if semantic_settings['Calculation'] == 'ProcessSearch':

            syntactic_settings.input.ams.PESExploration.Job = 'ProcessSearch'
            syntactic_settings.input.AMS.PESExploration.CalculateFragments = 'T'
            syntactic_settings.input.AMS.PESExploration.BindingSites.Calculate = 'T'
            syntactic_settings.input.AMS.PESExploration.DynamicSeedStates = 'T'

        syntactic_settings.input.ReaxFF.ForceField = \
            map_generic_setting(emmo.ForceFieldIdentifierString, root_cuds_object)

        syntactic_settings.input.ReaxFF.Charges.Solver = \
            map_generic_setting(emmo.Solver, root_cuds_object)

        syntactic_settings.input.AMS.Constraints.FixedRegion = \
            map_generic_setting(emmo.FixedRegion, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.NumExpeditions = \
            map_generic_setting(emmo.NumberOfExpeditions, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.NumExplorers = \
            map_generic_setting(emmo.NumberOfExplorers, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.SaddleSearch.MaxEnergy = \
            map_generic_setting(emmo.MaximumEnergy, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.StatesAlignment.ReferenceRegion = \
            map_generic_setting(emmo.ReferenceRegion, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.StructureComparison.CheckSymmetry = \
            map_generic_setting(emmo.CheckSymmetry, root_cuds_object)

        if map_generic_setting(emmo.CheckSymmetry, root_cuds_object) == 'T':

            # Hardcoded defaults
            syntactic_settings.input.AMS.PESExploration.StructureComparison.DistanceDifference = 0.1
            syntactic_settings.input.AMS.PESExploration.StructureComparison.EnergyDifference = 0.05
            syntactic_settings.input.AMS.PESExploration.StructureComparison.NeighborCutoff = 2.5

        syntactic_settings.input.AMS.PESExploration.RandomSeed = \
            map_generic_setting(emmo.RandomSeed, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.BindingSites.NeighborCutoff = \
            map_generic_setting(emmo.NeighborCutoff, root_cuds_object)

    elif semantic_settings['Calculation'] == 'BindingSites':
        previous = root_cuds_object.get(rel=emmo.hasSpatialNext.inverse).pop()

        syntactic_settings.input.AMS.task = "PESExploration"
        syntactic_settings.input.AMS.PESExploration.Job = 'BindingSites'
        syntactic_settings.input.AMS.PESExploration.BindingSites.Calculate = 'T'
        syntactic_settings.input.AMS.PESExploration.CalculateFragments = 'F'
        syntactic_settings.input.AMS.PESExploration.LoadEnergyLandscape.Path = \
            os.path.join(workdir, str(previous.uid))

        syntactic_settings.input.ReaxFF.ForceField = \
            map_generic_setting(emmo.ForceFieldIdentifierString, root_cuds_object)

        syntactic_settings.input.ReaxFF.Charges.Solver = \
            map_generic_setting(emmo.Solver, root_cuds_object)

        syntactic_settings.input.AMS.Constraints.FixedRegion = \
            map_generic_setting(emmo.FixedRegion, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.NumExpeditions = \
            map_generic_setting(emmo.NumberOfExpeditions, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.NumExplorers = \
            map_generic_setting(emmo.NumberOfExplorers, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.SaddleSearch.MaxEnergy = \
            map_generic_setting(emmo.MaximumEnergy, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.StatesAlignment.ReferenceRegion = \
            map_generic_setting(emmo.ReferenceRegion, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.StructureComparison.CheckSymmetry = \
            map_generic_setting(emmo.CheckSymmetry, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.RandomSeed = \
            map_generic_setting(emmo.RandomSeed, root_cuds_object)

        syntactic_settings.input.AMS.PESExploration.BindingSites.NeighborCutoff = \
            map_generic_setting(emmo.NeighborCutoff, root_cuds_object)

        # Hardcoded defaults
        if map_generic_setting(emmo.CheckSymmetry, root_cuds_object) == 'F':
            syntactic_settings.input.AMS.PESExploration.GenerateSymmetryImages = 'T'
            syntactic_settings.input.AMS.PESExploration.StructureComparison.\
                DistanceDifference = '0.1'
            syntactic_settings.input.AMS.PESExploration.StructureComparison.\
                EnergyDifference = '0.05'
            syntactic_settings.input.AMS.PESExploration.StructureComparison.\
                NeighborCutoff = '2.5'

    elif semantic_settings['Calculation'] == "StationaryPointCalculation":
        syntactic_settings.input.AMS.task = "TransitionStateSearch"

    elif semantic_settings['Calculation'] == "VibrationalFrequencyCalculation":
        syntactic_settings.input.AMS.task = "VibrationalAnalysis"

    # Map XC Functional:
    exceptions = ["PotentialEnergySurfaceExploration", "ProcessSearch",
                  "BindingSites", "LandscapeRefinement"]

    if semantic_settings['Calculation'] not in exceptions:

        semantic_settings['XCFunctional'] = map_xc_functional(
            root_cuds_object, semantic_settings['AccuracyLevel'], semantic_settings['Calculation'])

        # Map the ExchangeCorrelationFuntionals using AMS labels:
        LDA = ['Xonly', 'Xalpha', 'VWN', 'PW92']
        GGA = ['BP86', 'PW91', 'mPW', 'PBE', 'RPBE', 'revPBE', 'mPBE', 'PBEsol',
               'HTBS', 'BLYP', 'OLYP', 'OPBE', 'BEE', 'XLYP', 'SSB-D']
        HYBRID = ['B3LYP', 'B3LYP*', 'BHandHLYP', 'B1LYP', 'KMLYP', 'O3LYP',
                  'X3LYP', 'BHandH', 'B1PW91', 'mPW1PW', 'mPW1K', 'PBE0',
                  'OPBE0', 'S12H']

        if semantic_settings['XCFunctional'] in LDA:
            syntactic_settings.input.ADF.xc.LDA = semantic_settings['XCFunctional']

        elif semantic_settings['XCFunctional'] in GGA:
            syntactic_settings.input.ADF.xc.GGA = semantic_settings['XCFunctional']

        elif semantic_settings['XCFunctional'] in HYBRID:
            syntactic_settings.input.ADF.xc.HYBRID = semantic_settings['XCFunctional']

        else:
            raise_error(file=os.path.basename(__file__),
                        function=map_xc_functional.__name__,
                        type='NameError',
                        message="XC_functional label not implemented.")

        # Map Basis Set:
        semantic_settings['BasisSet'] = map_basis_set(root_cuds_object,
                                                      semantic_settings['AccuracyLevel'],
                                                      semantic_settings['Calculation'])
        syntactic_settings.input.ADF.basis.type = semantic_settings['BasisSet']
    return syntactic_settings


def map_accuracy_level(root_cuds_object: Cuds) -> str:
    """
    Map the GlobalAccuracy of the calculation, that will be used to set
    defaults.

    :param root_cuds_object: CUDS wrapper to be checked.

    :param engine: Engine provided when initializing the wrapper object.

    :param str: For the moment we just return a string.
    """

    search_accuracy_level = \
        search.find_cuds_objects_by_oclass(emmo.AccuracyLevel, root_cuds_object, emmo.hasInput)

    #  Raise error if two accuracies, for some reason, are found:
    if len(search_accuracy_level) > 1:
        raise_error(file=os.path.basename(__file__), function=map_accuracy_level.__name__,
                    type='NameError', message='More than one emmo.AccuracyLevel defined'
                    ' in the Wrapper object.')

    if search_accuracy_level:

        if search_accuracy_level.hasSymbolData.upper() in ['LOW', 'MEDIUM', 'HIGH']:
            return search_accuracy_level.hasSymbolData.upper()

    else:
        raise_warning(file=os.path.basename(__file__), function=map_accuracy_level.__name__,
                      message='No AccuracyLevel found, setting a default LOW accuracy')
    return 'LOW'


def map_calculation_type(root_cuds_object: Cuds) -> str:
    """Map the type of calculation of a CUDS Calculation object.

    :param root_cuds_object: CUDS Wrapper to be checked.

    :return str: The type of calculation. Sets WavefunctionOptimization as
                 default if any calculation is found in the wrapper.
    """

    search_calculation = \
        search.find_cuds_objects_by_oclass(emmo.Calculation, root_cuds_object, emmo.hasPart)

    search_simulation = \
        search.find_cuds_object(criterion=lambda x:
                                x.is_a(oclass=emmo.Simulation), root=root_cuds_object,
                                rel=emmo.hasPart, find_all=True, max_depth=1)

    #  Raise error if two calculations, for some reason, are found:
    if len(search_calculation) > 1 and len(search_simulation) == 0:
        raise_error(file=os.path.basename(__file__), function=map_calculation_type.__name__,
                    type='NameError', message='More than one emmo.calculation defined'
                                              ' in the Wrapper object.')

    if search_calculation:
        if 'WavefunctionOptimization' in str(search_calculation[0]):
            return 'WavefunctionOptimization'

        elif 'GeometryOptimization' in str(search_calculation[0]):
            return 'GeometryOptimization'

        elif 'PotentialEnergySurfaceScan' in str(search_calculation[0]):
            return 'PotentialEnergySurfaceScan'

        elif 'PotentialEnergySurfaceExploration' in str(search_calculation[0]):
            return 'PotentialEnergySurfaceExploration'

        elif 'ProcessSearch' in str(search_calculation[0]):
            return 'ProcessSearch'

        elif 'BindingSites' in str(search_calculation[0]):
            return 'BindingSites'

        elif 'StationaryPointCalculation' in str(search_calculation[0]):
            return 'StationaryPointCalculation'

        elif 'LandscapeRefinement' in str(search_calculation[0]):
            return 'LandscapeRefinement'

        elif 'VibrationalFrequencyCalculation' in str(search_calculation[0]):
            return 'VibrationalFrequencyCalculation'

        else:
            # If there a Calculation objet not defined above, set default:
            raise_warning(file=os.path.basename(__file__), function=map_calculation_type.__name__,
                          message='Calculation_type not implemented, setting'
                                  ' a default WavefunctionOptimization')
            return 'WavefunctionOptimization'

    return ''


def map_model_type(root_cuds_object: Cuds) -> str:
    """Map the model of calculation of a CUDS Calculation object.

    :param root_cuds_object: CUDS wrapper to be checked.

    :return str: The type of model.
    """

    search_calculation = \
        search.find_cuds_objects_by_oclass(emmo.Calculation, root_cuds_object, emmo.hasPart)

    search_simulation = \
        search.find_cuds_object(criterion=lambda x: x.is_a(oclass=emmo.Simulation),
                                root=root_cuds_object, rel=emmo.hasPart,
                                find_all=True, max_depth=1)

    #  Raise error if two calculations, for some reason, are found:
    if len(search_calculation) > 1 and len(search_simulation) == 0:
        raise_error(file=os.path.basename(__file__), function=map_calculation_type.__name__,
                    type='NameError',
                    message='More than one emmo.calculation defined'
                    ' in the Wrapper object.')

    if search_calculation:
        search_model = \
             search.find_cuds_objects_by_oclass(
                                emmo.DFTB,
                                search_calculation[0], emmo.hasInput)
        if len(search_model) == 1:
            return 'DFTB'
        elif len(search_model) > 1:
            raise_error(file=os.path.basename(__file__), function=map_PLAMSSettings.__name__,
                        type='NameError',
                        message='More than one Model CUDS object found in Calculation.')
        else:
            return ""
    else:
        return ""


def map_generic_setting(emmo_class: OntologyClass, root_cuds_object: Cuds, ) -> str:
    """Generic function to search for a parameter inside CUDS object.

    :param parameter: Element to search for in the CUDS object.

    :param root_cuds_object: CUDS wrapper to be checked.

    :return str: The equivalent string of the PLAMS setting.
    """

    search_emmo_class = \
        search.find_cuds_objects_by_oclass(
                           emmo_class,
                           root_cuds_object, emmo.hasPart)

    if len(search_emmo_class) > 1:
        raise_error(file=os.path.basename(__file__), function=map_generic_setting.__name__,
                    type='NameError', message='More than one equivalent PLAMSSetting defined'
                    ' in the Wrapper object.')

    if search_emmo_class:

        if search_emmo_class[0].is_a(emmo.FixedRegion):
            return search_emmo_class[0].hasSymbolData

        elif search_emmo_class[0].is_a(emmo.ForceFieldIdentifierString):

            if search_emmo_class[0].is_a(emmo.CHONSFPtClNi):
                return 'CHONSFPtClNi.ff'

            else:
                raise_error(file=os.path.basename(__file__), function=map_generic_setting.__name__,
                            type='NameError', message='ForceField is not implemented.')

        elif search_emmo_class[0].is_a(emmo.Solver):
            return (search_emmo_class[0].get(oclass=emmo.Symbol,
                                             rel=emmo.hasPart))[0].hasSymbolData

        elif search_emmo_class[0].is_a(emmo.NumberOfExpeditions):
            return search_emmo_class[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.NumberOfExplorers):
            return search_emmo_class[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.MaximumEnergy):
            return (search_emmo_class[0].get(oclass=emmo.Real,
                                             rel=emmo.hasPart))[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.ReferenceRegion):
            return search_emmo_class[0].hasSymbolData

        elif search_emmo_class[0].is_a(emmo.RandomSeed):
            return search_emmo_class[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.CheckSymmetry):
            return search_emmo_class[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.NeighborCutoff):
            return (search_emmo_class[0].get(oclass=emmo.Real,
                                             rel=emmo.hasPart))[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.ThermodynamicTemperature):
            return (search_emmo_class[0].get(oclass=emmo.Real,
                                             rel=emmo.hasPart))[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.Pressure):

            pressure = float(search_emmo_class[0].get(oclass=emmo.Real,
                                                      rel=emmo.hasPart)[0].hasNumericalData)
            pressure = pressure/100000

            return pressure

        elif search_emmo_class[0].is_a(emmo.MaximumTime):
            return (search_emmo_class[0].get(oclass=emmo.Real,
                                             rel=emmo.hasPart))[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.WallTime):
            return (search_emmo_class[0].get(oclass=emmo.Real,
                                             rel=emmo.hasPart))[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.MaximumSteps):
            return (search_emmo_class[0].get(oclass=emmo.Real,
                                             rel=emmo.hasPart))[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.Snapshots):
            return (search_emmo_class[0].get(oclass=emmo.Real,
                                             rel=emmo.hasSpatialPart))[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.ProcessStatistics):
            return (search_emmo_class[0].get(oclass=emmo.Real,
                                             rel=emmo.hasSpatialPart))[0].hasNumericalData

        elif search_emmo_class[0].is_a(emmo.SpeciesNumbers):
            return (search_emmo_class[0].get(oclass=emmo.Real,
                                             rel=emmo.hasSpatialPart))[0].hasNumericalData

    else:
        raise_error(file=os.path.basename(__file__),
                    function=map_generic_setting.__name__,
                    type='NameError',
                    message='PLAMSSetting is not implemented.')

    return


def map_xc_functional(root_cuds_object: Cuds, accuracy_level: str,
                      calculation_type: str) -> str:
    """
    Map ExchangeCorrelationFunctional in a Calculation CUDS object.

    If the container does not have the required object, it will set defaults
    according to accuracy_leven and calculation_type, searching in the dictionaries
    placed in:
    ./osp/defaults/dictionaries/

    If that happens, the XC_functional default will be added as CUDS object into the wrapper.

    :param root_cuds_object: = CUDS Calculation to be checked.

    :param accuracy_level: = str, accuracy level of the calculation.

    :param calculation_type: = str, type of calculation.

    :param return: String of the XC_functional.
    """

    available_xc_functionals = \
        emmo.ExchangeCorrelationFunctionalIdentifierString.direct_subclasses

    search_xc_functional = []

    for i in available_xc_functionals:
        xc_functional = search.find_cuds_object(criterion=lambda x: x.is_a(i),
                                                root=root_cuds_object, rel=emmo.hasPart,
                                                find_all=True)

        if xc_functional:
            search_xc_functional.extend(xc_functional)

    #  Raise error if two xc_functionals, for some reason, are found:
    if len(search_xc_functional) > 1:
        raise_error(file=os.path.basename(__file__), function=map_xc_functional.__name__,
                    type='NameError',
                    message='More than one emmo.ExchangeCorrelationFunctional'
                    ' defined in the Wrapper object.')

    if search_xc_functional:

        if 'B3LYP' in str(search_xc_functional[0]):
            return 'B3LYP'

        elif 'PBE' in str(search_xc_functional[0]):
            return 'PBE'

        elif 'VWN' or 'LDA' in str(search_xc_functional[0]):
            return 'VWN'

        else:
            raise_error(file=os.path.basename(__file__), function=map_xc_functional.__name__,
                        type='NameError', message='XC_functional selected is not implemented.')
    else:
        default_setting = AMS_default_setting(root_cuds_object, accuracy_level,
                                              calculation_type, 'ExchangeCorrelationFunctional')
    return default_setting


def map_force_field(root_cuds_object: Cuds, engine: str) -> str:
    """
    Map ForceField in a wrapper CUDS object.

    There is no default as this option will be triggered only if a PESExploration
    is inside the CUDS object.

    :param root_cuds_object: = CUDS wrapper to be checked.

    :param engine: Engine provided when initializing the wrapper object.

    :param str: String of the ForceField
    """

    search_force_field = \
        search.find_cuds_objects_by_oclass(
                           emmo.ForceFieldIdentifierString,
                           root_cuds_object, emmo.hasPart)

    if search_force_field:
        if 'AMBER' in str(search_force_field[0]):
            return 'AMBER'

        elif 'ANAKINME' in str(search_force_field[0]):
            return 'ANAKINME'

        elif 'CHONSFPtClNi' in str(search_force_field[0]):
            return 'CHONSFPtClNi'

        elif 'OPLS' in str(search_force_field[0]):
            return 'OPLS'

        elif 'OPLSAA' in str(search_force_field[0]):
            return 'OPLSAA'

        else:
            raise_error(file=os.path.basename(__file__), function=map_force_field.__name__,
                        type='NameError', message='Force Field selected is not implemented.')
    return


def map_basis_set(root_cuds_object: Cuds, global_accuracy: str,
                  calculation_type: str) -> str:
    """
    Map BasisSet in a Calculation CUDS object.

    If the container does not have the required object, it will set defaults
    according to dictionaries in:
    ./osp/defaults/dictionaries/

    The BasisSet default will be added as CUDS object into the Calculation object.

    :param root_cuds_object: CUDS wrapper to be checked.

    :param global_accuracy: Global accuracy of the workflow.

    :param calculation_type: Calculation type set already in the wrapper.

    :param str: String of the BasisSet.
    """

    available_xc_basis_sets = \
        emmo.BasisSetIdentifierString.direct_subclasses

    search_basis_set = []

    for i in available_xc_basis_sets:
        basis_set = search.find_cuds_object(criterion=lambda x: x.is_a(i), root=root_cuds_object,
                                            rel=emmo.hasPart, find_all=True)
        if basis_set:
            search_basis_set.extend(basis_set)

    # Raise error if two basis_set, for some reason, are found:
    if len(search_basis_set) > 1:
        raise_error(file=os.path.basename(__file__), function=map_basis_set.__name__,
                    type='NameError', message='More than one emmo.BasisSet defined'
                    ' in the Wrapper object.')

    if search_basis_set:
        if 'SZ' in str(search_basis_set[0]):
            return 'SZ'
        elif 'DZP' in str(search_basis_set[0]):
            return 'DZP'
        elif 'DZ' in str(search_basis_set[0]):
            return 'DZ'
        elif 'TZP' in str(search_basis_set[0]):
            return 'TZP'
        elif 'TZ2P' in str(search_basis_set[0]):
            return 'TZ2P'
        else:
            raise_error(file=os.path.basename(__file__), function=map_basis_set.__name__,
                        type='NameError', message='BasisSet selected is not implemented.')
    else:
        default_setting = AMS_default_setting(root_cuds_object,
                                              global_accuracy,
                                              calculation_type, 'BasisSet')
    return default_setting


def map_PyZacrosSettings(root_cuds_object: Cuds) -> pz.Settings:
    """
    Map PyZacrosSettings from CUDS objects to engine.

    :param root_cuds_object: Calculation CUDS object.

    :return pz.Settings: The mapped settings.
    """

    syntactic_settings = pz.Settings()

    syntactic_settings.Random_Seed = map_generic_setting(emmo.RandomSeed, root_cuds_object)

    syntactic_settings.Temperature = map_generic_setting(emmo.ThermodynamicTemperature,
                                                         root_cuds_object)

    syntactic_settings.Pressure = map_generic_setting(emmo.Pressure, root_cuds_object)

    # Search molar fractions:
    # PyZacros list molar fractions as settings this is why the search is done here.
    search_molar_fraction = search.find_cuds_objects_by_oclass(emmo.AmountFraction,
                                                               root_cuds_object, emmo.hasInput)
    if search_molar_fraction:

        molar_fraction_dict = {'molar_fraction': {}}
        for molar_fraction in search_molar_fraction:
            atom_label = molar_fraction.get(oclass=emmo.ChemicalElement,
                                            rel=emmo.hasProperty)[0].hasSymbolData
            molar_fraction_value = molar_fraction.get(oclass=emmo.Real,
                                                      rel=emmo.hasQuantityValue)[0].hasNumericalData
            molar_fraction_dict['molar_fraction'][atom_label] = float(molar_fraction_value)

        molar_fraction_settings = PlamsSettings(molar_fraction_dict)
        syntactic_settings = syntactic_settings + molar_fraction_settings

    else:
        # look deeper into the species attached, especially useful when input is coming from
        # a BaseModel
        search_gas_species = search.find_cuds_objects_by_oclass(emmo.GasSpecies,
                                                                root_cuds_object, emmo.hasInput)

        molar_fraction_dict = {'molar_fraction': {}}
        for gas_species in search_gas_species:
            molar_fraction = gas_species.get(oclass=emmo.AmountFraction, rel=emmo.hasProperty)
            molar_fraction_value = molar_fraction[0].get(
                oclass=emmo.Real, rel=emmo.hasQuantityValue)[0].hasNumericalData
            atom_label = gas_species.get(oclass=emmo.ChemicalElement,
                                         rel=emmo.hasProperty)[0].hasSymbolData
            molar_fraction_dict['molar_fraction'][atom_label] = float(molar_fraction_value)

        molar_fraction_settings = PlamsSettings(molar_fraction_dict)
        syntactic_settings = syntactic_settings + molar_fraction_settings

    # Search snapshots:
    search_snapshots = search.find_cuds_objects_by_oclass(emmo.Snapshots,
                                                          root_cuds_object, emmo.hasInput)

    if search_snapshots:
        snapshots = float(map_generic_setting(emmo.Snapshots, root_cuds_object))
        search_snapshots[0].hasSymbolData.replace("on ", "")
        syntactic_settings.snapshots = (search_snapshots[0].hasSymbolData.replace("on ", ""),
                                        snapshots)

    # Search process_statistics:
    search_process_statistics = \
        search.find_cuds_objects_by_oclass(emmo.ProcessStatistics, root_cuds_object, emmo.hasInput)

    if search_process_statistics:
        process_statistics = float(map_generic_setting(emmo.ProcessStatistics, root_cuds_object))
        syntactic_settings.process_statistics = (
            search_process_statistics[0].hasSymbolData.replace("on ", ""), process_statistics)

    # Search species_numbers:
    search_species_numbers = \
        search.find_cuds_objects_by_oclass(emmo.SpeciesNumbers, root_cuds_object, emmo.hasInput)

    if search_species_numbers:
        species_numbers = float(map_generic_setting(emmo.SpeciesNumbers, root_cuds_object))
        syntactic_settings.species_numbers = (
            search_species_numbers[0].hasSymbolData.replace("on ", ""), species_numbers)

    # Search max_steps:
    search_max_steps = search.find_cuds_objects_by_oclass(emmo.MaximumSteps, root_cuds_object,
                                                          emmo.hasInput)

    if search_max_steps:
        max_time = map_generic_setting(emmo.MaximumSteps, root_cuds_object)

        if max_time == "inf":
            syntactic_settings.Max_Steps = "infinity"

        else:
            syntactic_settings.Max_Steps = max_time

    # Search max_time:
    search_max_time = \
        search.find_cuds_objects_by_oclass(emmo.MaximumTime, root_cuds_object, emmo.hasInput)

    if search_max_time:
        max_time = map_generic_setting(emmo.MaximumTime, root_cuds_object)
        if max_time == "inf":
            syntactic_settings.Max_Time = "infinity"
        else:
            syntactic_settings.Max_Time = max_time

    # Search wall_time:
    search_wall_time = search.find_cuds_objects_by_oclass(emmo.WallTime, root_cuds_object,
                                                          emmo.hasInput)

    if search_wall_time:
        syntactic_settings.Wall_Time = map_generic_setting(emmo.WallTime, root_cuds_object)

    return syntactic_settings

def find_reaction_path(root_cuds_object: Cuds) -> List[Cuds]:
    """
    Find reaction pathways in root cuds object with specific order

    :param root_cuds_object: Calculation CUDS object.

    :return List[Cuds]:
    """
    landscape = search.find_cuds_objects_by_oclass(emmo.EnergyLandscape, root_cuds_object, emmo.hasPart)
    if len(landscape) > 1:
        raise_error(file=os.path.basename(__file__),
                    function=find_reaction_path.__name__,
                    type='ValueError',
                    message='More than one energy landscape found.')
    elif not landscape:
        raise_error(file=os.path.basename(__file__),
                    function=find_reaction_path.__name__,
                    type='ValueError',
                    message='No energy landscape found.')
    else:
        first = landscape.pop().get(oclass=emmo.ReactionPathway, rel=emmo.hasSpatialFirst)
        if len(first) > 1:
            raise_error(file=os.path.basename(__file__),
                        function=find_reaction_path.__name__,
                        type='ValueError',
                        message='More than one first spatial reaction pathways found.')
        elif not first:
            raise_error(file=os.path.basename(__file__),
                        function=find_reaction_path.__name__,
                        type='ValueError',
                        message='First spatial pathway in energy landscape not found.')
        else:
            paths = []
            current = first
            while current:
                current = current.pop()
                paths.append(current)
                current = current.get(rel=emmo.hasSpatialNext)
            return paths

def map_PyZacrosMechanism(root_cuds_object: Cuds) -> pz.Mechanism:
    """
    Map ontology-defined Mechanism to pz.Mechanism() object.

    :param root_cuds_object: Calculation CUDS object.

    :return PyZacros Mechanism:
    """

    search_reversible_reaction = \
        search.find_cuds_object(criterion=lambda x:
                                x.is_a(oclass=emmo.ReversibleReactionEquation),
                                root=root_cuds_object, rel=emmo.hasPart,
                                find_all=True, max_depth=1)

    if search_reversible_reaction:

        reaction_list = list()

        for index, equation in enumerate(search_reversible_reaction):

            # Label:
            if "mechanism_input" in str(equation.uid):
                reaction_label = equation.uid.split("step.")[1]
            else:
                reaction_label = "Reaction" + str(1+index)

            # Find initial Species:
            search_reactants = \
                search.find_cuds_object(criterion=lambda x:
                                        x.is_a(oclass=emmo.ChemicalReactionEquationReactant),
                                        root=equation, rel=emmo.hasPart,
                                        find_all=True, max_depth=1)

            # Find final Species:
            search_products = \
                search.find_cuds_object(criterion=lambda x:
                                        x.is_a(oclass=emmo.ChemicalReactionEquationProduct),
                                        root=equation, rel=emmo.hasPart,
                                        find_all=True, max_depth=1)

            # Find Gas Reactant Species:
            search_gas_reactants = \
                search.find_cuds_object(criterion=lambda x:
                                        x.is_a(oclass=emmo.GasReactantSpecies),
                                        root=equation, rel=emmo.hasPart,
                                        find_all=True, max_depth=1)

            # Find Gas Product Species:
            search_gas_products = \
                search.find_cuds_object(criterion=lambda x:
                                        x.is_a(oclass=emmo.GasProductSpecies),
                                        root=equation, rel=emmo.hasPart,
                                        find_all=True, max_depth=1)

            # Find Sites:
            search_site_types = \
                search.find_cuds_object(criterion=lambda x: x.is_a(oclass=crystallography.Site),
                                        root=equation, rel=emmo.hasPart, find_all=True, max_depth=1)

            if search_site_types:
                site_types = []

                for site in search_site_types:
                    site_types.append(site.get(oclass=emmo.ChemicalElement,
                                               rel=emmo.hasSpatialDirectPart)[0].hasSymbolData)
            else:
                site_types = None

            # Find neighboring:
            search_neighbors = \
                search.find_cuds_object(criterion=lambda x: x.is_a(oclass=emmo.Neighboring),
                                        root=equation, rel=emmo.hasPart, find_all=True, max_depth=1)

            if search_neighbors:
                neighbor_list = []

                # only implemented in the case of two neighbors:
                # -1 to match Zacros to Python-based indexing
                neighbor_tuple = (
                    int(search_neighbors[0].get(
                        oclass=emmo.Real, rel=emmo.hasSpatialFirst)[0].hasNumericalData)-1,
                    int(search_neighbors[0].get(oclass=emmo.Real,
                                                rel=emmo.hasSpatialNext)[0].hasNumericalData)-1)
                neighbor_list.append(neighbor_tuple)

            else:
                neighbor_list = None

            # Find pre_expon:
            search_pre_expon = \
                search.find_cuds_object(criterion=lambda x:
                                        x.is_a(oclass=emmo.ArrheniusCoefficient), root=equation,
                                        rel=emmo.hasPart, find_all=True, max_depth=1)

            if search_pre_expon:
                pre_expon = float(search_pre_expon[0].get(oclass=emmo.Real,
                                  rel=emmo.hasSpatialPart)[0].hasNumericalData)

            # Find pre_expon:
            search_pe_ratio = \
                search.find_cuds_object(criterion=lambda x:
                                        x.is_a(oclass=emmo.Constant), root=equation,
                                        rel=emmo.hasPart, find_all=True, max_depth=1)

            if search_pe_ratio:
                pe_ratio = float(search_pe_ratio[0].get(
                    oclass=emmo.Real, rel=emmo.hasSpatialPart)[0].hasNumericalData)
            else:
                pe_ratio = 0.0

            # Find activ_eng:
            search_activ_eng = \
                search.find_cuds_object(criterion=lambda x: x.is_a(oclass=emmo.ActivationEnergy),
                                        root=equation, rel=emmo.hasPart,
                                        find_all=True, max_depth=1)

            if search_activ_eng:
                activ_eng = float(search_activ_eng[0].get(
                    oclass=emmo.Real, rel=emmo.hasSpatialPart)[0].hasNumericalData)

        # Map reactants:
            if search_reactants:
                reactant_list = []
                for reactant in search_reactants:
                    reactant_list.append(reactant.get(oclass=emmo.ChemicalElement,
                                         rel=emmo.hasSpatialPart)[0].hasSymbolData)

            if reactant_list:
                reactant_object_list = []
                reactant_object_list = [pz.Species(symbol=reactant) for reactant in reactant_list]

            if search_gas_reactants:
                gas_reactant_list = []
                gas_reactant_formation_energy_list = []

                for reactant in search_gas_reactants:

                    gas_reactant_list.append(reactant.get(oclass=emmo.ChemicalElement,
                                             rel=emmo.hasSpatialPart)[0].hasSymbolData)
                    search_formation_energy = \
                        search.find_cuds_object(criterion=lambda x:
                                                x.is_a(oclass=emmo.FormationEnergy), root=reactant,
                                                rel=emmo.hasPart, find_all=True, max_depth=1)

                    reactant_formation_energy_value = search_formation_energy[0].get(
                        oclass=emmo.Real, rel=emmo.hasPart)[0].hasNumericalData

                    gas_reactant_formation_energy_list.append(
                        float(reactant_formation_energy_value))

                gas_reactant_object_list = []

                for index, reactant in enumerate(gas_reactant_list):
                    gas_reactant_object_list.append(pz.Species(
                        symbol=reactant, gas_energy=gas_reactant_formation_energy_list[index]))

                reactant_object_list = reactant_object_list + gas_reactant_object_list

        # Map products:
            if search_products:

                product_list = []
                for product in search_products:
                    product_list.append(product.get(oclass=emmo.ChemicalElement,
                                        rel=emmo.hasSpatialPart)[0].hasSymbolData)
            if product_list:

                product_object_list = []
                product_object_list = [pz.Species(symbol=product) for product in product_list]

            if search_gas_products:
                gas_product_list = []
                gas_product_formation_energy_list = []

                for product in search_gas_products:

                    gas_product_list.append(
                        product.get(oclass=emmo.ChemicalElement,
                                    rel=emmo.hasSpatialPart)[0].hasSymbolData)
                    search_formation_energy = \
                        search.find_cuds_object(criterion=lambda x:
                                                x.is_a(oclass=emmo.FormationEnergy),
                                                root=product, rel=emmo.hasPart,
                                                find_all=True, max_depth=1)
                    product_formation_energy_value = search_formation_energy[0].get(
                        oclass=emmo.Real, rel=emmo.hasPart)[0].hasNumericalData
                    gas_product_formation_energy_list.append(float(product_formation_energy_value))

                gas_product_object_list = []

                for index, product in enumerate(gas_product_list):
                    gas_product_object_list.append(
                        pz.Species(symbol=product,
                                   gas_energy=gas_product_formation_energy_list[index]))

                product_object_list = product_object_list + gas_product_object_list

            reaction = pz.ElementaryReaction(
                initial=reactant_object_list, final=product_object_list,
                reversible=True, pre_expon=pre_expon, pe_ratio=pe_ratio,
                neighboring=neighbor_list, activation_energy=activ_eng,
                label=reaction_label, site_types=site_types)
            neighbor_list = []
            reaction_list.append(reaction)

        syntactic_mechanism = pz.Mechanism(reaction_list)

    else:
        raise_error(file=os.path.basename(__file__),
                    function=map_PyZacrosMechanism.__name__,
                    type='NameError',
                    message='No ReversibleReactionEquation found in ReactionMechanism.')

    return syntactic_mechanism


def map_PyZacrosClusterExpansion(semantic_cluster_list: list) -> pz.ClusterExpansion:
    """
    Map ontology-defined ClusterExpansion list to pz.ClusterExpansion() object.

    :param cluster_list: List of CLusterExpansion CUDS objects

    :return PyZacros ClusterExpansion:
    """

    syntactic_cluster_list = []

    for index, cluster in enumerate(semantic_cluster_list):

        species_list = []
        denticity_list = []
        entity_number_list = []
        neighbor_list = []

        # Label:
        if "energetics_input.name." in str(cluster.uid):
            cluster_label = cluster.uid.split("energetics_input.name.")[1]
        else:
            cluster_label = "Cluster" + str(1+index)

        # Find Sites:
        search_site_types = \
            search.find_cuds_object(criterion=lambda x: x.is_a(oclass=crystallography.Site),
                                    root=cluster, rel=emmo.hasPart, find_all=True, max_depth=1)

        if search_site_types:

            site_types_list = []

            for site in search_site_types:
                site_types_list.append(site.get(oclass=emmo.ChemicalElement,
                                       rel=emmo.hasSpatialDirectPart)[0].hasSymbolData)
        else:
            site_types_list = None

        # Find neighboring:
        search_neighbors = \
            search.find_cuds_object(criterion=lambda x: x.is_a(oclass=emmo.Neighboring),
                                    root=cluster, rel=emmo.hasSpatialPart, find_all=True,
                                    max_depth=1)

        if search_neighbors:

            if len(site_types_list) == 2:
                # special case of two neighbors:
                # -1 to match Zacros to Python-based indexing
                neighbor_tuple = (
                    int(search_neighbors[0].get(oclass=emmo.Real,
                        rel=emmo.hasSpatialFirst)[0].hasNumericalData)-1,
                    int(search_neighbors[0].get(oclass=emmo.Real,
                        rel=emmo.hasSpatialNext)[0].hasNumericalData)-1)
                neighbor_list.append(neighbor_tuple)
        else:
            neighbor_list = None

        # Finds species
        search_lattice_states = \
            search.find_cuds_object(criterion=lambda x: x.is_a(oclass=emmo.LatticeState),
                                    root=cluster, rel=emmo.hasPart, find_all=True, max_depth=1)

        for lattice_state in search_lattice_states:
            search_species = \
                search.find_cuds_object(criterion=lambda x: x.is_a(oclass=emmo.AdsorbedSpecies),
                                        root=lattice_state, rel=emmo.hasPart, find_all=True,
                                        max_depth=1)

            for species in search_species:
                species_list.append(species.get(oclass=emmo.ChemicalElement,
                                    rel=emmo.hasSpatialDirectPart)[0].hasSymbolData)

                entity_number_list.append(
                    int(species.get(oclass=emmo.Real,
                                    rel=emmo.hasSpatialDirectPart)[0].hasNumericalData)-1)

                search_denticity = \
                    search.find_cuds_object(criterion=lambda x: x.is_a(oclass=emmo.Denticity),
                                            root=species, rel=emmo.hasPart,
                                            find_all=True, max_depth=1)

                denticity_list.append(int(search_denticity[0].get(oclass=emmo.Real,
                                      rel=emmo.hasSpatialDirectPart)[0].hasNumericalData))

        # Fix denticity for bidentate species:
        if len(species_list) == 2:
            if species_list[0] == species_list[1]:
                if entity_number_list[0] == entity_number_list[1]:
                    denticity_list[0] = max(denticity_list[0], denticity_list[1])
                    denticity_list[1] = max(denticity_list[0], denticity_list[1])

        # Find cluster_energy:
        search_cluster_energy = \
            search.find_cuds_object(criterion=lambda x: x.is_a(oclass=emmo.ClusterEnergy),
                                    root=cluster, rel=emmo.hasPart, find_all=True, max_depth=1)

        if search_cluster_energy:
            cluster_energy = float(search_cluster_energy[0].get(oclass=emmo.Real,
                                   rel=emmo.hasSpatialPart)[0].hasNumericalData)

        # Generation of the syntactic_objects
        species_object_list = []
        for index, species in enumerate(species_list):
            species_object_list.append(pz.Species(symbol=species,
                                                  denticity=denticity_list[index]))

        # Map individual clusters:
        syntactic_cluster = pz.Cluster(
                            species=species_object_list, site_types=site_types_list,
                            entity_number=entity_number_list, neighboring=neighbor_list,
                            energy=cluster_energy, label=cluster_label)

        syntactic_cluster_list.append(syntactic_cluster)

    return pz.ClusterExpansion(syntactic_cluster_list)


def map_energies(root_cuds_object: Cuds) -> list:
    """
    Gather a list of three TotalElectronicEnergy values from a cuba.Wrapper
    object.
    """

    calculations = []
    electronic_energies = []
    electronic_energy_values = []
    search_simulation = \
        search.find_cuds_object(criterion=lambda x: x.is_a(oclass=emmo.Simulation),
                                root=root_cuds_object, rel=emmo.hasPart,
                                find_all=True, max_depth=1)

    if not len(search_simulation):
        return
    
    # First
    next_calc = search_simulation[0].get(rel=emmo.hasSpatialFirst)

    while next_calc:
        current = next_calc.pop()
        calculations.append(current)
        next_calc = current.get(rel=emmo.hasSpatialNext)

    for calculation in calculations:

        electronic_energy = \
            search.find_cuds_objects_by_oclass(oclass=emmo.TotalElectronicEnergy,
                                               root=calculation, rel=emmo.hasOutput)
        electronic_energies += electronic_energy

    for te in electronic_energies:

        real_number = te.get(oclass=emmo.Real)
        electronic_energy_values.append(real_number[0].hasNumericalData)

    return electronic_energy_values


def map_results(engine, root_cuds_object: Cuds) -> str:
    """
    Add main results to the Cuds (Wrapper) object.

    :engine: main engine object, e.g. Multijob for AMS calculations.
    :param root_cuds_object: CUDS wrapper to add results to.
    :param wrapper_settings: SimamsSettings object.
    :return str: file system path to tarball-
   """

    if isinstance(engine, MultiJob):

        # in case of several children, outputs will be added in the order of dependencies.
        search_simulation = search.find_cuds_object(
            criterion=lambda x: x.is_a(oclass=emmo.Simulation) or x.is_a(oclass=emmo.Workflow),
            root=root_cuds_object,
            rel=cuba.relationship,
            find_all=True,
            max_depth=1
        )

        if not search_simulation:

            # Attach results to Calculation CUDS object using hasOutput rel.
            search_calculation = \
                    search.find_cuds_objects_by_oclass(
                                       emmo.Calculation,
                                       root_cuds_object, rel=cuba.relationship)
            tarball = map_tarball(engine, search_calculation[0])

            if map_calculation_type(root_cuds_object) == "WavefunctionOptimization":
                add_AMSenergy_to_object(engine.children[0], search_calculation[0])

            elif map_calculation_type(root_cuds_object) == "GeometryOptimization":

                optimized_geometry = read_molecule(
                    engine.children[0].results.get_molecule(section="Molecule"),
                    GO_optimized=True, read_from_object=True)
                search_calculation[0].add(optimized_geometry, rel=emmo.hasOutput)
                add_AMSenergy_to_object(engine.children[0], search_calculation[0])

            elif map_calculation_type(root_cuds_object) == "ProcessSearch":

                # Print some results as output
                print(engine.children[0].results.get_energy_landscape())

                # Attach Mechanism and Cluster_Expansions to Wrapper object as Output:
                loader_ads = pz.RKFLoader(engine.children[0].results)
                loader_ads.replace_site_types(['A', 'B', 'C'], ['fcc', 'br', 'hcp'])

                # 1. Mechanism
                mechanism_output = read_mechanism(str(loader_ads.mechanism), read_from_file=False)
                search_calculation[0].add(mechanism_output, rel=emmo.hasOutput)

                # 2. Cluster
                cluster_output = read_cluster_expansion(str(loader_ads.clusterExpansion),
                                                        read_from_file=False)
                for cluster in cluster_output:
                    search_calculation[0].add(cluster, rel=emmo.hasOutput)

                print(loader_ads.clusterExpansion)
                print(loader_ads.mechanism)

            elif map_calculation_type(root_cuds_object) == "BindingSites":

                # Attach UnitCell() with file path to Wrapper object as Output:
                loader_bs = pz.RKFLoader(engine.children[0].results)
                loader_bs.replace_site_types(['A', 'B', 'C'], ['fcc', 'br', 'hcp'])
                loader_bs.lattice.set_repeat_cell((10, 10))
                loader_bs.lattice.plot()
                
                file = os.path.join(engine.path, "lattice_input.dat")
                with open(file, "w+") as lattice_file:
                    lattice_file.write(str(loader_bs.lattice))
                uuid = get_upload(file)
                lattice_output = crystallography.UnitCell(uid=UUID(uuid))
                search_calculation[0].add(lattice_output, rel=emmo.hasOutput)

            elif map_calculation_type(root_cuds_object) == "LandscapeRefinement":
                add_AMSLandscape_to_object(engine.children[0], search_calculation[0])

            # TODO elif map_calculation_type(root_cuds_object) == "XXXX":

        else:

            simulation = search_simulation[0]

            # First
            next_calc = simulation.get(rel=emmo.hasSpatialFirst)
            i = 0

            while next_calc:
                current = next_calc.pop()
                if current.is_a(emmo.WavefunctionOptimization) \
                        or current.is_a(emmo.GeometryOptimization):
                    add_AMSenergy_to_object(engine.children[i], current)

                elif current.is_a(emmo.ProcessSearch):
                    # Print some results as output
                    print(engine.children[i].results.get_energy_landscape())

                    # Attach Mechanism and Cluster_Expansions to Wrapper object as Output:
                    loader_ads = pz.RKFLoader(engine.children[i].results)
                    loader_ads.replace_site_types(['A', 'B', 'C'], ['fcc', 'br', 'hcp'])

                    # 1. Mechanism
                    mechanism_output = read_mechanism(str(loader_ads.mechanism), read_from_file=False)
                    current.add(mechanism_output, rel=emmo.hasOutput)
                    if simulation.is_a(emmo.Simulation):
                        simulation.add(mechanism_output, rel=emmo.hasOutput)

                    # 2. Cluster
                    cluster_output = read_cluster_expansion(str(loader_ads.clusterExpansion),
                                                            read_from_file=False)
                    for cluster in cluster_output:
                        current.add(cluster, rel=emmo.hasOutput)
                        if simulation.is_a(emmo.Simulation):
                            simulation.add(cluster, rel=emmo.hasOutput)

                    print(loader_ads.clusterExpansion)
                    print(loader_ads.mechanism)
                elif current.is_a(emmo.BindingSites):
                    # Attach UnitCell() with file path to Wrapper object as Output:
                    loader_bs = pz.RKFLoader(engine.children[i].results)
                    loader_bs.replace_site_types(['A', 'B', 'C'], ['fcc', 'br', 'hcp'])
                    loader_bs.lattice.set_repeat_cell((10, 10))
                    loader_bs.lattice.plot()

                    file = os.path.join(engine.path, "lattice_input.dat")
                    with open(file, "w+") as lattice_file:
                        lattice_file.write(str(loader_bs.lattice))
                    uuid = get_upload(file)
                    lattice_output = crystallography.UnitCell(uid=UUID(uuid))
                    current.add(lattice_output, rel=emmo.hasOutput)
                    if simulation.is_a(emmo.Simulation):
                        simulation.add(lattice_output, rel=emmo.hasOutput)

                i += 1
                next_calc = current.get(oclass=emmo.AtomisticCalculation, rel=emmo.hasSpatialNext) \
                    or current.get(oclass=emmo.PostProcessing, rel=emmo.hasSpatialNext)

            if simulation.is_a(emmo.Simulation):
                tarball = map_tarball(engine, simulation)
            elif simulation.is_a(emmo.Workflow):
                tarball = map_tarball(engine, current)

    elif isinstance(engine, pz.ZacrosJob):
        search_calculation = search.find_cuds_objects_by_oclass(
                                   emmo.MesoscopicCalculation, root_cuds_object,
                                   rel=cuba.relationship)
        tarball = map_tarball(engine, search_calculation[0])

    else:
        raise_error(file=os.path.basename(__file__), function=map_results.__name__,
                    type='NameError',
                    message='Mapping of the results is not implemented for this engine.')
    return tarball



def add_AMSenergy_to_object(job: AMSJob, root_cuds_object: Cuds):
    """
    Add Total Energy result to the Cuds (Calculation) object.

    :job: Multijob object.

    :param root_cuds_object: CUDS Calculation to add results to.
   """

    energy = emmo.TotalElectronicEnergy()
    energy_value = emmo.Real(hasNumericalData=job.results.get_energy(unit='eV'))
    energy_unit = emmo.ElectronVolt(hasSymbolData='eV')
    energy.add(energy_value, rel=emmo.hasQuantityValue)
    energy.add(energy_unit, rel=emmo.hasReferenceUnit)
    root_cuds_object.add(energy, rel=emmo.hasOutput)
#    print("Total Energy value:", job.results.get_energy(unit='eV'), "eV")

    return


def add_AMSLandscape_to_object(job: AMSJob, root_cuds_object: Cuds):
    """
    Add EnergyLandscape to the Cuds (Calculation) object.

    :job: Multijob object.

    :param root_cuds_object: CUDS Calculation to add results to.
   """

    syntactic_energy_landscape = job.results.get_energy_landscape()
    semantic_energy_landscape = emmo.EnergyLandscape()

    for state in syntactic_energy_landscape:

        semantic_path = emmo.ReactionPathway()
        semantic_molecule = read_molecule(filename=state.molecule, read_from_object=True)

        # Add TotalEnergy
        energy = emmo.TotalElectronicEnergy()
        energy_value = emmo.Real(
                hasNumericalData=state.energy)
        energy_unit = emmo.HartreeEnergy(hasSymbolData='Ha')
        energy.add(energy_value, rel=emmo.hasQuantityValue)
        energy.add(energy_unit, rel=emmo.hasReferenceUnit)
        semantic_molecule.add(energy, rel=emmo.hasProperty)

        # Add Charge
        if state.molecule.properties.charge:

            semantic_molecule_charge = emmo.ElectricCharge()
            semantic_molecule_integer = emmo.Integer(
                hasNumericalData=state.molecule.properties.charge)
            semantic_molecule_charge.add(semantic_molecule_integer, rel=emmo.hasQuantityValue)
            semantic_molecule.add(semantic_molecule_charge, rel=emmo.hasProperty)

        if state.isTS:

            # ## REACTANTS ##
            semantic_reactant = emmo.ChemicalReactionEquationReactant()
            semantic_reactant_molecule = read_molecule(filename=state.reactants.molecule,
                                                       read_from_object=True)

            # Add TotalEnergy
            energy = emmo.TotalElectronicEnergy()
            energy_value = emmo.Real(hasNumericalData=state.reactants.energy)
            energy.add(energy_value, rel=emmo.hasQuantityValue)
            semantic_reactant_molecule.add(energy, rel=emmo.hasProperty)

            # Add Charge
            if state.reactants.molecule.properties.charge:
                semantic_reactant_molecule_charge = emmo.ElectricCharge()
                semantic_reactant_molecule_integer = emmo.Integer(
                    hasNumericalData=state.reactants.molecule.properties.charge)
                semantic_reactant_molecule_charge.add(
                    semantic_reactant_molecule_integer, rel=emmo.hasQuantityValue)
                semantic_reactant_molecule.add(
                    semantic_reactant_molecule_charge, rel=emmo.hasProperty)

            semantic_reactant.add(semantic_reactant_molecule, rel=emmo.hasSpatialDirectPart)
            semantic_path.add(semantic_reactant, rel=emmo.hasPart)

            # ## PRODUCTS ##
            semantic_product = emmo.ChemicalReactionEquationProduct()
            semantic_product_molecule = read_molecule(
                filename=state.products.molecule, read_from_object=True)

            # Add TotalEnergy
            energy = emmo.TotalElectronicEnergy()
            energy_value = emmo.Real(hasNumericalData=state.products.energy)
            energy.add(energy_value, rel=emmo.hasQuantityValue)
            semantic_product_molecule.add(energy, rel=emmo.hasProperty)

            # Add Charge
            if state.products.molecule.properties.charge:
                semantic_product_molecule_charge = emmo.ElectricCharge()
                semantic_product_molecule_integer = emmo.Integer(
                    hasNumericalData=state.products.molecule.properties.charge)
                semantic_product_molecule_charge.add(
                    semantic_product_molecule_integer, rel=emmo.hasQuantityValue)
                semantic_product_molecule.add(semantic_product_molecule_charge,
                                              rel=emmo.hasProperty)

            semantic_product.add(semantic_product_molecule, rel=emmo.hasSpatialDirectPart)
            semantic_path.add(semantic_product, rel=emmo.hasPart)

            # ## TS ##
            semantic_TS = emmo.TransitionStateGeometry()
            semantic_TS.add(semantic_molecule, rel=emmo.hasSpatialDirectPart)
            semantic_path.add(semantic_TS, rel=emmo.hasPart)

        else:
            semantic_path.add(semantic_molecule, rel=emmo.hasPart)

        semantic_energy_landscape.add(semantic_path, rel=emmo.hasPart)

    root_cuds_object.add(semantic_energy_landscape, rel=emmo.hasOutput)
#    pretty_print(semantic_energy_landscape)

    return


def map_tarball(engine, root_cuds_object, path=None) -> str:
    """
    Creates and links tarball of the results to root_cuds_object Calculation:

    :engine: Job PLAMS /Multijob PLAMS object for AMS calculations, str: for Zacros calculation

    :param root_cuds_object: CUDS Calculation to add results to.

    :return str: file system path to tarball
   """

    tar = tempfile.NamedTemporaryFile().name
    shutil.make_archive(tar, "tar", engine.path)
    tar_file = f"{tar}.tar"

    print(f"Job output dumped to {tar_file}")

    iri = arcp_random(tar_file)
    string = emmo.String(iri=iri)
    root_cuds_object.add(string, rel=emmo.hasOutput)
    return tar_file
