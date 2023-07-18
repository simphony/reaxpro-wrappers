from pathlib import Path
from scm.plams import Molecule as PlamsMolecule
from scm.plams import Atom as PlamsAtom
from osp.tools.io_functions import molecule_from_PlamsMolecule, \
                                   molecule_from_xyz, lattice_from_xyz, read_mechanism, \
                                   read_cluster_expansion
from osp.tools.mapping_functions import map_molecule, map_lattice, \
                                        map_PyZacrosMechanism, map_PyZacrosClusterExpansion


EXAMPLES = Path(__file__).absolute().parents[1] / 'tests/test_files'


def test_molecule_from_xyz():
    """ Test function reading a .xyz formatted file."""
    filename = EXAMPLES / 'CO_ads+Pt111.xyz'
    dict_to_test = {}
    dict_molecule_xyz = {
      'atom_symbol':
      ['C', 'O', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt',
       'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt',
       'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt',
       'Pt', 'Pt', 'Pt', 'Pt', 'Pt'],
      'atom_coordinates':
      [['5.08349688', '2.95549819', '9.52500529'],
       ['5.08420252', '2.95701878', '10.74456106'],
       ['9.24070201', '6.95527267', '1.31036083'],
       ['9.24070201', '2.15427277', '1.31036083'],
       ['10.6266313', '4.55477272', '1.31036083'],
       ['12.01256059', '6.95527267', '1.31036083'],
       ['3.69698485', '2.15427277', '8.1'],
       ['5.08291414', '4.55477272', '8.1'],
       ['6.46884343', '6.95527267', '8.1'],
       ['6.46884343', '2.15427277', '8.1'],
       ['7.85477272', '4.55477272', '8.1'],
       ['9.24070201', '6.95527267', '8.1'],
       ['9.24070201', '2.15427277', '8.1'],
       ['10.6266313', '4.55477272', '8.1'],
       ['12.01256059', '6.95527267', '8.1'],
       ['5.08291414', '2.95443942', '5.83678694'],
       ['6.46884343', '5.35493937', '5.83678694'],
       ['7.85477272', '7.75543932', '5.83678694'],
       ['7.85477272', '2.95443942', '5.83678694'],
       ['9.24070201', '5.35493937', '5.83678694'],
       ['10.6266313', '7.75543932', '5.83678694'],
       ['10.6266313', '2.95443942', '5.83678694'],
       ['12.01256059', '5.35493937', '5.83678694'],
       ['13.39848987', '7.75543932', '5.83678694'],
       ['2.31105556', '1.35410612', '3.57357389'],
       ['3.69698485', '3.75460607', '3.57357389'],
       ['5.08291414', '6.15510602', '3.57357389'],
       ['5.08291414', '1.35410612', '3.57357389'],
       ['6.46884343', '3.75460607', '3.57357389'],
       ['7.85477272', '6.15510602', '3.57357389'],
       ['7.85477272', '1.35410612', '3.57357389'],
       ['9.24070201', '3.75460607', '3.57357389'],
       ['10.6266313', '6.15510602', '3.57357389'],
       ['3.69698485', '2.15427277', '1.31036083'],
       ['5.08291414', '4.55477272', '1.31036083'],
       ['6.46884343', '6.95527267', '1.31036083'],
       ['6.46884343', '2.15427277', '1.31036083'],
       ['7.85477272', '4.55477272', '1.31036083']],
      'region': ['adsorbate', 'adsorbate', 'surface', 'surface',
                 'surface', 'surface', 'surface', 'surface',
                 'surface', 'surface', 'surface', 'surface',
                 'surface', 'surface', 'surface', 'surface',
                 'surface', 'surface', 'surface', 'surface',
                 'surface', 'surface', 'surface', 'surface',
                 'surface', 'surface', 'surface', 'surface',
                 'surface', 'surface', 'surface', 'surface',
                 'surface', 'surface', 'surface', 'surface',
                 'surface', 'surface']}

    test_molecule = molecule_from_xyz(filename, TS=False)
    dict_to_test = map_molecule(test_molecule)

    assert dict_to_test == dict_molecule_xyz


def test_molecule_from_PlamsMolecule():
    """ Test the mapping of a PlamsMolecule."""
    syntactic_molecule = PlamsMolecule()
    syntactic_molecule.add_atom(PlamsAtom(symbol='O', coords=(0, 0, 0)))
    syntactic_molecule.add_atom(PlamsAtom(symbol='H', coords=(1, 0, 0)))
    syntactic_molecule.add_atom(PlamsAtom(symbol='H', coords=(0, 1, 0)))

    dict_to_test = map_molecule(molecule_from_PlamsMolecule(syntactic_molecule, GO_optimized=False))

    dict_PlamsMolecule = {'atom_symbol':
                          ['O', 'H', 'H'],
                          'atom_coordinates':
                          [['0.0', '0.0', '0.0'],
                              ['1.0', '0.0', '0.0'],
                              ['0.0', '1.0', '0.0']]}

    assert dict_to_test == dict_PlamsMolecule


def test_lattice_from_xyz():
    """ Test lattice by reading a .xyz formatted file."""
    filename = EXAMPLES / 'CO_ads+Pt111.xyz'
    dict_to_test = {}

    dict_lattice = {'vector_coordinates': [(8.31557575, 0.0, 0.0), (4.15778787, 7.20149984, 0.0)]}

    test_molecule = lattice_from_xyz(filename)
    dict_to_test = map_lattice(test_molecule)

    assert dict_to_test == dict_lattice


def test_PZ_mechanism() -> str:
    """Test reading and mapping pyZacros Mechanisms from Zacros formatted files."""
    filename = EXAMPLES / 'mechanism_input.dat'
    str_to_test = str(map_PyZacrosMechanism(read_mechanism(filename)))

    str_mechanism = \
"""mechanism

reversible_step Reaction1
  gas_reacs_prods O2 -1
  sites 2
  neighboring 1-2
  initial
    1 * 1
    2 * 1
  final
    1 O* 1
    2 O* 1
  site_types brg brg
  pre_expon  7.98000e+07
  pe_ratio  9.43100e-09
  activ_eng  0.00000e+00
end_reversible_step

reversible_step Reaction2
  gas_reacs_prods CO -1
  sites 1
  initial
    1 * 1
  final
    1 CO* 1
  site_types brg
  pre_expon  4.26500e+07
  pe_ratio  6.56300e-09
  activ_eng  0.00000e+00
end_reversible_step

reversible_step Reaction3
  gas_reacs_prods CO2 1
  sites 2
  neighboring 1-2
  initial
    1 CO* 1
    2 O* 1
  final
    1 * 1
    2 * 1
  site_types brg brg
  pre_expon  2.78600e+12
  pe_ratio  3.23100e+07
  activ_eng  5.20000e-01
end_reversible_step

end_mechanism"""
    assert str_to_test == str_mechanism


def test_PZ_cluster():
    """Test reading and mapping pyZacros Clusters from Zacros formatted files."""
    filename = EXAMPLES / 'energetics_input.dat'
    str_to_test = str(map_PyZacrosClusterExpansion(read_cluster_expansion(filename)))
    str_cluster = \
"""energetics

cluster Cluster1
  sites 1
  lattice_state
    1 CO* 1
  site_types brg
  graph_multiplicity 1
  cluster_eng -2.36000e+00
end_cluster

cluster Cluster2
  sites 1
  lattice_state
    1 CO* 1
  site_types hol
  graph_multiplicity 1
  cluster_eng -1.85000e+00
end_cluster

cluster Cluster3
  sites 1
  lattice_state
    1 O* 1
  site_types brg
  graph_multiplicity 1
  cluster_eng -1.51000e+00
end_cluster

cluster Cluster4
  sites 2
  neighboring 1-2
  lattice_state
    1 CO* 1
    2 O* 1
  site_types hol hol
  graph_multiplicity 1
  cluster_eng  5.00000e-02
end_cluster

cluster Cluster5
  sites 2
  neighboring 1-2
  lattice_state
    1 CO2* 1
    1 CO2* 2
  site_types brg hol
  graph_multiplicity 1
  cluster_eng -3.42000e+00
end_cluster

end_energetics"""

    assert str_to_test == str_cluster
