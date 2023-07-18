"""AMS/EON calculation using ReaxPro ontology."""

from osp.core.namespaces import emmo, cuba
from osp.wrappers.simams.simams_session import SimamsSession
from osp.tools.io_functions import read_molecule, read_lattice
#from osp.core.utils import pretty_print

calculation = emmo.ProcessSearch()

molecule = read_molecule('./XYZ/CO_ads+Pt111.xyz')
lattice = read_lattice('./XYZ/CO_ads+Pt111.xyz')

forcefield = emmo.CHONSFPtClNi()

solver = emmo.Solver()
type_of_solver = emmo.Symbol(hasSymbolData='Direct') 
solver.add(type_of_solver, rel=emmo.hasPart)

fixed_region = emmo.FixedRegion(hasSymbolData='surface') 
num_expeditions = emmo.NumberOfExpeditions(hasNumericalData=30)
num_explorers = emmo.NumberOfExplorers(hasNumericalData=4)

max_energy = emmo.MaximumEnergy()
energy_value = emmo.Real(hasNumericalData='2.0')
energy_unit = emmo.ElectronVolt(hasSymbolData='eV')
max_energy.add(energy_value, rel=emmo.hasQuantityValue)
max_energy.add(energy_unit, rel=emmo.hasReferenceUnit)

ref_region = emmo.ReferenceRegion(hasSymbolData='surface') 
random_seed = emmo.RandomSeed(hasNumericalData='100')

max_distance = emmo.NeighborCutoff()
distance_value = emmo.Real(hasNumericalData='3.8')
distance_unit = emmo.Ångström(hasSymbolData='Å')
max_distance.add(distance_value, rel=emmo.hasQuantityValue)
max_distance.add(distance_unit, rel=emmo.hasReferenceUnit)


symmetry_check = emmo.CheckSymmetry(hasNumericalData = 'T')

calculation.add(molecule, lattice, forcefield, solver, fixed_region, num_expeditions,
                num_explorers, max_energy, ref_region, random_seed, max_distance, symmetry_check, 
                rel=emmo.hasInput)

with SimamsSession() as sess:
    reaxpro_wrapper = cuba.Wrapper(session=sess)
    reaxpro_wrapper.add(calculation, rel=emmo.hasPart)
    reaxpro_wrapper.session.run()