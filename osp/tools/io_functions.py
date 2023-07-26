"""Support functions to write and read data."""
import os
from enum import Enum
from pathlib import Path
import yaml

from scm.plams import Molecule as PlamsMolecule
from scm.plams.interfaces.adfsuite.ams import AMSResults

from osp.core.cuds import Cuds
from osp.core.namespaces import crystallography, emmo
# from osp.core.utils import pretty_print


class GeometryType(str, Enum):

    XYZ = "xyz"
    MOL = "mol"
    PDB = "pdb"


def read_molecule(filename: str, TS: bool = False, GO_optimized: bool = False,
                  read_from_object: bool = False, geometry_type: GeometryType = None):
    """
    Read geometry from an external file or from PLAMS object.

    Supported formats: xyz, mol and pdb.

    :param filename: Name of the file to be read.

    :param TS: wether the output CUDS object should be a emmo.TransitionStateGeometry or not.

    :param GO_optimized: in case the output be a emmo.OptimizedMolecularGeometry.

    :param read_from_object: Optionally, emmo.MolecularGeometry is retrieved from a
                             PLAMSMolecule object.

    :param return: emmo.MolecularGeometry or emmo.TransitionStateGeometry or
                   emmo.OptimizedMolecularGeometry CUDS object.
    """

    if not read_from_object:
        if not geometry_type:
            geometry_type = assign_geom_type(filename)

        if geometry_type == GeometryType.XYZ:
            return molecule_from_xyz(filename, TS)

        elif geometry_type == GeometryType.MOL:
            return molecule_from_mol(filename, TS)

        elif geometry_type == GeometryType.PDB:
            raise_error(
                file=os.path.basename(__file__),
                function=read_molecule.__name__,
                type="NotImplementedError",
                message="pdb format not implemented yet.",
            )
    else:
        # Read from PLAMS.molecule format.
        return molecule_from_PlamsMolecule(filename, GO_optimized)


def read_lattice(filename: str) -> crystallography.UnitCell:
    """
    Function to read lattice vectors from different file formats.

    Supported formats: xyz, mol and pdb.

    :param filename: Name of the file to be read.

    :param return: crystallography.UnitCell object with added LatticeVectors.
    """

    geometry_type = assign_geom_type(filename)

    if geometry_type == "xyz":
        return lattice_from_xyz(filename)

    elif geometry_type == "mol":
        raise_error(
            file=os.path.basename(__file__),
            function=read_molecule.__name__,
            type="NotImplementedError",
            message="mol format not implemented yet.",
        )

    elif geometry_type == "pdb":
        raise_error(
            file=os.path.basename(__file__),
            function=read_molecule.__name__,
            type="NotImplementedError",
            message="pdb format not implemented yet.",
        )
    return


def assign_geom_type(filename: str) -> str:
    """
    Assign a geometry type format to a entry filename file.

    :param filename: name containing the path and name to file.

    :param return: String with the format of the geometry file.
    """

    path = Path(filename)

    suffix = path.suffix.lower()

    if suffix not in (".xyz", ".mol", ".pdb"):
        raise NameError(f"Unknown geometry format: {suffix}.")

    return suffix.strip(".")


def assign_chemical_name(root_cuds_object: Cuds, labels: list) -> Cuds:
    """
    Assign a ChemicalName to a Molecular Geometry.

    :param root_cuds_object: emmo.Molecule

    :param labels: list with atom labels.

    :param return: same CUDS object with the ChemicalName added.
    """
    counts = []
    labels_no_redundant = []
    label = str()

    for lb in labels:
        counts.append(labels.count(lb))

    for index, lb in enumerate(labels):
        if counts[index] == 1:
            labels[index] = lb
        else:
            labels[index] = lb + str(counts[index])

    # Eliminate redundancy on labels:
    labels_no_redundant = sorted(list(dict.fromkeys(labels)))

    # Get the final molecular label:
    for lb in labels_no_redundant:
        if lb[:3] != "VEC":
            label = label + lb

    chemical_name = emmo.ChemicalName(hasSymbolData=label)
    root_cuds_object.add(chemical_name, rel=emmo.hasSign)

    return


def molecule_from_xyz(filename: str, TS: bool):
    """
    Create emmo.MolecularGeometry or emmo.TransitionStateGeometry from .xyz file.

    It adds single CUDS atom containers with CUDS atom_labels and CUDS positions.

    :param filename: name containing the path and name to file.

    :param TS: option to return emmo.TransitionStateGeometry instead of emmo.MolecularGeometry.

    :param return: emmo.MolecularGeometry containing atoms, labels and position objects.
    """

    if TS:
        transition_state = emmo.TransitionStateGeometry()
    else:
        molecule = emmo.MolecularGeometry()

    labels = []
    coordinates = []
    regions = []
    xyz = open(filename)
    n_atoms = int(xyz.readline())
    _ = xyz.readline()

    # Reading the .xyz file:
    for line in xyz:

        params = line.split()
        label, x, y, z = params[0:4]
        try:
            if params[4]:
                region = params[4].strip()
        except IndexError:
            region = "x=NIL"

        labels.append(label)
        coordinates.append([float(x), float(y), float(z)])
        regions.append(region.split("=")[1])
    xyz.close()
    print("molecule_from_xyz: External", filename, "file has been read.")

    # Creating atom containers and adding them to the root_cuds_object
    # (i.e. MolecularGeometry)
    atoms = []
    for iatom in range(n_atoms):
        atom = emmo.AtomEntity()
        atom_symbol = emmo.ChemicalElement(hasSymbolData=labels[iatom])

        position_vec = emmo.PositionVector()
        coordinate_x = emmo.Real(hasNumericalData=coordinates[iatom][0])
        coordinate_y = emmo.Real(hasNumericalData=coordinates[iatom][1])
        coordinate_z = emmo.Real(hasNumericalData=coordinates[iatom][2])
        position_vec.add(coordinate_x, rel=emmo.hasSpatialFirst)
        position_vec.add(coordinate_y, rel=emmo.hasSpatialNext)
        position_vec.add(coordinate_z, rel=emmo.hasSpatialLast)

        if regions[iatom] == "NIL":
            atom.add(position_vec, atom_symbol, rel=emmo.hasPart)
        else:
            region = emmo.Region(hasSymbolData=regions[iatom].lower())
            atom.add(position_vec, atom_symbol, region, rel=emmo.hasPart)
        if iatom == 0:
            rel = emmo.hasSpatialFirst
        elif iatom == n_atoms:
            rel = emmo.hasSpatialLast
        else:
            rel = emmo.hasSpatialDirectPart
            atoms[iatom-1].add(atom, rel=emmo.hasSpatialNext)
        atoms.append(atom)
        if TS:
            transition_state.add(atom, rel=rel)
        else:
            molecule.add(atom, rel=rel)
    if TS:
        assign_chemical_name(transition_state, labels)
        return transition_state
    else:
        assign_chemical_name(molecule, labels)
        return molecule


def molecule_from_mol(filename: str, TS: bool = False):
    """
    Create emmo.MolecularGeometry from .mol file.

    It adds single CUDS atom containers with CUDS atom_labels and CUDS positions.

    :param filename: name containing the path and name to file.

    :param TS: option to return emmo.TransitionStateGeometry instead of emmo.MolecularGeometry.

    :param return: emmo.MolecularGeometry or emmo.TransitionStateGeometry containing atoms,
                   labels and position objects.

    #TODO, info on the Bonding structure is not read.
    """

    if TS:
        transition_state = emmo.TransitionStateGeometry()
    else:
        molecule = emmo.MolecularGeometry()

    labels = []
    coordinates = []
    with open(filename, "r") as file:
        mol = file.readlines()
        n_atoms = int(mol[3].split()[0])

        # Reading atoms from .mol file:
        for line in mol[4: 4 + n_atoms]:

            params = line.split()
            x, y, z, label = params[0:4]

            labels.append(label)
            coordinates.append([float(x), float(y), float(z)])

    print("molecule_from_mol External", filename, "file has been read.")

    # Creating atom containers and adding them to the root_cuds_object
    # (i.e. MolecularGeometry)
    for iatom in range(n_atoms):
        atom = emmo.AtomEntity()
        atom_symbol = emmo.ChemicalElement(hasSymbolData=labels[iatom])

        position_vec = emmo.PositionVector()
        coordinate_x = emmo.Real(hasNumericalData=coordinates[iatom][0])
        coordinate_y = emmo.Real(hasNumericalData=coordinates[iatom][1])
        coordinate_z = emmo.Real(hasNumericalData=coordinates[iatom][2])
        position_vec.add(coordinate_x, rel=emmo.hasSpatialFirst)
        position_vec.add(coordinate_y, rel=emmo.hasSpatialNext)
        position_vec.add(coordinate_z, rel=emmo.hasSpatialLast)

        atom.add(position_vec, atom_symbol, rel=emmo.hasPart)
        if TS:
            transition_state.add(atom, rel=emmo.hasSpatialPart)
        else:
            molecule.add(atom, rel=emmo.hasSpatialPart)
    if TS:
        assign_chemical_name(transition_state, labels)
        return transition_state
    else:
        assign_chemical_name(molecule, labels)
        return molecule


def molecule_from_PlamsMolecule(syntactic_molecule: PlamsMolecule, GO_optimized: bool):
    """
    Create emmo.MolecularGeometry from a string from a PlamsMolecule.

    It adds single CUDS atom containers with CUDS atom_labels and CUDS positions.

    :param filename: bulk string.

    :param GO_optimized: in case the output be a emmo.OptimizedMolecularGeometry.

    :param return: emmo.MolecularGeometry or emmo.OptimizedMolecularGeometry containing atoms,
                   labels and position objects.
    """
    labels = []

    if GO_optimized:
        optimized_geometry = emmo.OptimizedMolecularGeometry()
    else:
        molecule = emmo.MolecularGeometry()

    # Reading the .xyz file:
    for atom in syntactic_molecule:

        labels.append(atom.symbol)

        semantic_atom = emmo.AtomEntity()

        atom_symbol = emmo.ChemicalElement(hasSymbolData=atom.symbol)
        position_vec = emmo.PositionVector()
        coordinate_x = emmo.Real(hasNumericalData=atom.coords[0])
        coordinate_y = emmo.Real(hasNumericalData=atom.coords[1])
        coordinate_z = emmo.Real(hasNumericalData=atom.coords[2])
        position_vec.add(coordinate_x, rel=emmo.hasSpatialFirst)
        position_vec.add(coordinate_y, rel=emmo.hasSpatialNext)
        position_vec.add(coordinate_z, rel=emmo.hasSpatialLast)
        semantic_atom.add(position_vec, atom_symbol, rel=emmo.hasPart)

        if GO_optimized:
            optimized_geometry.add(semantic_atom, rel=emmo.hasSpatialPart)
        else:
            molecule.add(semantic_atom, rel=emmo.hasSpatialPart)

    if GO_optimized:
        assign_chemical_name(optimized_geometry, labels)
        return optimized_geometry
    else:
        assign_chemical_name(molecule, labels)
        return molecule


def lattice_from_xyz(filename: str) -> crystallography.UnitCell:
    """
    Create a list of emmo.LatticeVectors from .xyz file.

    :param filename: name containing the path and name to file.

    :param return: List of emmo.LatticeVectors.
    """

    lattice = crystallography.UnitCell()
    labels = []
    coordinates = []
    regions = []
    xyz = open(filename, "r")
    n_atoms = int(xyz.readline())
    _ = xyz.readline()
    n_lines = 2  # to account for the second .xyz blank line.

    # Reading the .xyz file:
    for line in xyz:

        params = line.split()
        label, x, y, z = params[0:4]
        try:
            if params[4]:
                region = params[4].strip()
        except IndexError:
            region = "x=NIL"

        labels.append(label)
        coordinates.append([float(x), float(y), float(z)])
        regions.append(region.split("=")[1])
        n_lines += 1
    xyz.close()
    print("lattice_from_xyz: External", filename, "file has been read.")

    n_vectors = n_lines - n_atoms - 2  # n_atoms line + blank line

    # Raise error if no vectors
    if n_vectors == 0:
        raise_error(
            file=os.path.basename(__file__),
            function=lattice_from_xyz.__name__,
            type="NameError",
            message="LatticeVectors not found in .xyz file.",
        )

    for ivector in range(n_vectors):
        lattice_vec = crystallography.LatticeVector()
        coordinate_x = emmo.Real(
            hasNumericalData=coordinates[len(coordinates) - 1 - ivector][0]
        )
        coordinate_y = emmo.Real(
            hasNumericalData=coordinates[len(coordinates) - 1 - ivector][1]
        )
        coordinate_z = emmo.Real(
            hasNumericalData=coordinates[len(coordinates) - 1 - ivector][2]
        )
        lattice_vec.add(coordinate_x, rel=emmo.hasSpatialFirst)
        lattice_vec.add(coordinate_y, rel=emmo.hasSpatialNext)
        lattice_vec.add(coordinate_z, rel=emmo.hasSpatialLast)
        if ivector == 0:
            lattice.add(lattice_vec, rel=emmo.hasSpatialNext)
        if ivector == 1:
            lattice.add(lattice_vec, rel=emmo.hasSpatialFirst)
        if ivector == 2:
            lattice.add(lattice_vec, rel=emmo.hasSpatialLast)
    return lattice


def read_mechanism(
    filename: str, read_from_file: bool = True
) -> emmo.ChemicalReactionMechanism:
    """
    Create emmo.ChemicalReactionMechanism from a Zacros-formatted .dat file.

    :param filename: name containing the path and name to file.

    :param return: emmo.ChemicalReactionMechanism with all the semantic information.
    """

    mechanism = emmo.ChemicalReactionMechanism()

    if read_from_file is True:
        dat = open(filename, "r")
    elif read_from_file is False:
        dat = (
            filename.splitlines()
        )  # For cases where the file names contains the full string.

    initial = False
    final = False

    for count, line in enumerate(dat):
        params = line.split()
        if params:
            if params[0] == "reversible_step":
                reaction = emmo.ReversibleReactionEquation()

            if params[0] == "gas_reacs_prods":
                filedir = os.path.dirname(os.path.abspath(__file__))
                dictionary = os.path.join(
                    filedir, "../dictionaries/energies/gas_energies.yaml"
                )
                with open(dictionary, "r") as f:
                    energy_list = yaml.load(f, Loader=yaml.FullLoader)
                    formation_energy_value_from_dic = energy_list["Gas_energies"][
                        params[1]
                    ]
                if float(params[2]) < 0:
                    gas_reactant = emmo.GasReactantSpecies()
                    atom_symbol = emmo.ChemicalElement(hasSymbolData=params[1])
                    formation_energy = emmo.FormationEnergy()
                    formation_energy_value = emmo.Real(
                        hasNumericalData=formation_energy_value_from_dic
                    )
                    formation_energy.add(
                        formation_energy_value, rel=emmo.hasQuantityValue
                    )
                    coefficient = emmo.StoichiometricCoefficient()
                    coefficient_value = emmo.Real(hasNumericalData=params[2])
                    coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                    gas_reactant.add(
                        atom_symbol,
                        coefficient,
                        formation_energy,
                        rel=emmo.hasSpatialDirectPart
                    )
                    reaction.add(gas_reactant, rel=emmo.hasSpatialPart)
                if float(params[2]) > 0:
                    gas_reactant = emmo.GasProductSpecies()
                    atom_symbol = emmo.ChemicalElement(hasSymbolData=params[1])
                    formation_energy = emmo.FormationEnergy()
                    formation_energy_value = emmo.Real(
                        hasNumericalData=formation_energy_value_from_dic
                    )
                    formation_energy.add(
                        formation_energy_value, rel=emmo.hasQuantityValue
                    )
                    coefficient = emmo.StoichiometricCoefficient()
                    coefficient_value = emmo.Real(hasNumericalData=params[2])
                    coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                    gas_reactant.add(
                        atom_symbol,
                        coefficient,
                        formation_energy,
                        rel=emmo.hasSpatialDirectPart,
                    )
                    reaction.add(gas_reactant, rel=emmo.hasSpatialPart)

            if params[0] == "sites":
                num_of_sites = int(params[1])

            if params[0] == "neighboring":
                if num_of_sites == 2:
                    neighbors = emmo.Neighboring()
                    for i in range(num_of_sites):
                        neighbor_site = emmo.Real(
                            hasNumericalData=params[1].split("-")[i]
                        )
                        if i == 0:
                            neighbors.add(neighbor_site, rel=emmo.hasSpatialFirst)
                        else:
                            neighbors.add(neighbor_site, rel=emmo.hasSpatialNext)
                    reaction.add(neighbors, rel=emmo.hasSpatialPart)

            if params[0] == "initial":
                record = count
                initial = True

            if initial:
                if count == record:
                    continue
                elif count <= record + num_of_sites:
                    reactant = emmo.ChemicalReactionEquationReactant()
                    atom_symbol = emmo.ChemicalElement(hasSymbolData=params[1])
                    coefficient = emmo.StoichiometricCoefficient()
                    coefficient_value = emmo.Real(hasNumericalData=params[2])
                    coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                    reactant.add(
                        atom_symbol, coefficient, rel=emmo.hasSpatialDirectPart
                    )
                    reaction.add(reactant, rel=emmo.hasSpatialPart)
                else:
                    initial = False
                    #

            if params[0] == "final":
                record = count
                final = True

            if final:
                if count == record:
                    continue
                elif count <= record + num_of_sites:
                    reactant = emmo.ChemicalReactionEquationProduct()
                    atom_symbol = emmo.ChemicalElement(hasSymbolData=params[1])
                    coefficient = emmo.StoichiometricCoefficient()
                    coefficient_value = emmo.Real(hasNumericalData=params[2])
                    coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                    reactant.add(
                        atom_symbol, coefficient, rel=emmo.hasSpatialDirectPart
                    )
                    reaction.add(reactant, rel=emmo.hasSpatialPart)
                else:
                    final = False
                    #

            if params[0] == "site_types":
                site_type = crystallography.Site()
                atom_symbol = emmo.ChemicalElement(hasSymbolData=params[1])
                site_type.add(atom_symbol, rel=emmo.hasSpatialDirectPart)
                reaction.add(site_type, rel=emmo.hasSpatialPart)
                if num_of_sites == 2:
                    site_type = crystallography.Site()
                    atom_symbol = emmo.ChemicalElement(hasSymbolData=params[2])
                    site_type.add(atom_symbol, rel=emmo.hasSpatialDirectPart)
                    reaction.add(site_type, rel=emmo.hasSpatialPart)

            if params[0] == "pre_expon":
                coefficient = emmo.ArrheniusCoefficient()
                coefficient_value = emmo.Real(hasNumericalData=params[1])
                coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                reaction.add(coefficient, rel=emmo.hasSpatialPart)

            if params[0] == "pe_ratio":
                coefficient = emmo.Constant()
                coefficient_value = emmo.Real(hasNumericalData=params[1])
                coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                reaction.add(coefficient, rel=emmo.hasSpatialPart)

            if params[0] == "activ_eng":
                coefficient = emmo.ActivationEnergy()
                coefficient_value = emmo.Real(hasNumericalData=params[1])
                coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                reaction.add(coefficient, rel=emmo.hasSpatialPart)

            if params[0] == "end_reversible_step":
                mechanism.add(reaction, rel=emmo.hasPart)

    #    pretty_print(mechanism)

    if read_from_file is True:
        dat.close()

    print("read_mechanism: External", filename, "file has been read.")

    return mechanism


def read_cluster_expansion(filename: str, read_from_file: bool = True) -> list:
    """
    Create a [emmo.ClusterExpansion] list from a Zacros-formatted .dat file.

    :param filename: name containing the path and name to file.

    :param return: list of emmo.ClusterExpansion
    """

    if read_from_file is True:
        dat = open(filename, "r")
    elif read_from_file is False:
        dat = (
            filename.splitlines()
        )  # For cases where the file names contains the full string.

    cluster_list = []
    read_lattice_state = False

    for count, line in enumerate(dat):
        params = line.split()
        if params:

            if params[0] == "cluster":
                cluster = emmo.ClusterExpansion()

            if params[0] == "sites":
                num_of_sites = int(params[1])

            if params[0] == "neighboring":
                if num_of_sites == 2:
                    neighbors = emmo.Neighboring()
                    for i in range(num_of_sites):
                        neighbor_site = emmo.Real(
                            hasNumericalData=params[1].split("-")[i]
                        )
                        if i == 0:
                            neighbors.add(neighbor_site, rel=emmo.hasSpatialFirst)
                        else:
                            neighbors.add(neighbor_site, rel=emmo.hasSpatialNext)
                    cluster.add(neighbors, rel=emmo.hasSpatialPart)

            if params[0] == "lattice_state":
                record = count
                read_lattice_state = True

            if read_lattice_state:
                if count == record:
                    lattice_state = emmo.LatticeState()
                    continue
                elif count <= record + num_of_sites:
                    adsorbed_species = emmo.AdsorbedSpecies()
                    species_number = emmo.Real(hasNumericalData=params[0])
                    atom_symbol = emmo.ChemicalElement(hasSymbolData=params[1])

                    denticity_number = emmo.Real(hasNumericalData=params[2])
                    denticity = emmo.Denticity()
                    denticity.add(denticity_number, rel=emmo.hasQuantityValue)

                    adsorbed_species.add(
                        species_number,
                        atom_symbol,
                        denticity,
                        rel=emmo.hasSpatialDirectPart,
                    )

                    lattice_state.add(adsorbed_species, rel=emmo.hasSpatialDirectPart)
                else:
                    cluster.add(lattice_state, rel=emmo.hasPart)
                    read_lattice_state = False

            if params[0] == "site_types":
                site_type = crystallography.Site()
                atom_symbol = emmo.ChemicalElement(hasSymbolData=params[1])
                site_type.add(atom_symbol, rel=emmo.hasSpatialDirectPart)
                cluster.add(site_type, rel=emmo.hasSpatialPart)
                if num_of_sites == 2:
                    site_type = crystallography.Site()
                    atom_symbol = emmo.ChemicalElement(hasSymbolData=params[2])
                    site_type.add(atom_symbol, rel=emmo.hasSpatialDirectPart)
                    cluster.add(site_type, rel=emmo.hasSpatialPart)

            if params[0] == "cluster_eng":
                cluster_energy = emmo.ClusterEnergy()
                cluster_energy_value = emmo.Real(hasNumericalData=params[1])
                cluster_energy.add(cluster_energy_value, rel=emmo.hasQuantityValue)
                cluster.add(cluster_energy, rel=emmo.hasSpatialPart)

            if params[0] == "end_cluster":
                cluster_list.append(cluster)

    #    for i in cluster_list:
    #        pretty_print(i)

    if read_from_file is True:
        dat.close()

    print("read_cluster_expansion: External", filename, "file has been read.")

    return cluster_list


def print_results(plams_results: AMSResults, result: str):
    """
    Print results of a PLAMS/AMS/ADF calculation.

    :param plams_results: PLAMS result object to be analyzed.

    :param result: String to select what information to print.
    """

    if result == "total_energy":
        total_energy = plams_results.get_energy(unit="eV")
        print("Total electronic energy:", total_energy, "eV")

    if result == "energy_per_go_step":
        print("Total energy per geometry optimization step:")
        nEntries = plams_results.readrkf("History", "nEntries")
        for k in range(nEntries):
            s = "Energy({})".format(k + 1)
            x = plams_results.readrkf("History", s)
            print("Optimization cycle", k + 1, ":", x, "Hartree")
    return


def raise_error(file: str, function: str, type: str, message: str):
    """
    Raise errors according to arguments:

    :param file: Python file where the error occurs.

    :param function: Name of the function calling the error.

    :param type: Error type.

    :param message: String describing the error.
    """

    nl = "\n"

    msg = (
        f"{nl}### ERROR ### {function} function/method"
        f" in {file}:"
        f'{nl}              "{message}"'
    )
    if type == "NameError":
        raise NameError(msg)
    elif type == "ValueError":
        raise ValueError(msg)
    elif type == "NotImplementedError":
        raise NameError(message)

def raise_warning(file: str, function: str, message: str, message2=None):
    """
    Raise warning according to arguments:

    :param file: Python file where the error occurs.

    :param function: Name of the function calling the error.

    :param message: String describing the error.
    """

    nl = "\n"
    if not message2:
        msg = (
            f"{nl}### WARNING! ### {function} function/method"
            f" in {file}:"
            f'{nl}              "{message}"'
        )
    else:
        msg = (
            f"{nl}### WARNING! ### {function} function/method"
            f" in {file}:"
            f'{nl}              "{message}:", {message2}'
        )
    print(msg)
    return
