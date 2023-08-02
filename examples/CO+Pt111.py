from osp.core.namespaces import emmo, cuba, crystallography
from osp.wrappers.simams.simams_session import SimamsSession
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession 
from osp.tools.io_functions import read_molecule, read_lattice
from osp.core.utils import simple_search as search
from osp.core.utils import pretty_print, Cuds2dot, export_cuds

# PES Exploration 

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

max_distance = emmo.NeighborCutoff()
distance_value = emmo.Real(hasNumericalData='3.8')
distance_unit = emmo.Ångström(hasSymbolData='Å')
max_distance.add(distance_value, rel=emmo.hasQuantityValue)
max_distance.add(distance_unit, rel=emmo.hasReferenceUnit)

ref_region = emmo.ReferenceRegion(hasSymbolData='surface') 
random_seed = emmo.RandomSeed(hasNumericalData='100')

symmetry_check = emmo.CheckSymmetry(hasNumericalData = 'T')
calculation.add(molecule, lattice, solver, fixed_region, num_expeditions,
                forcefield, num_explorers, max_energy, max_distance, ref_region,
                random_seed, symmetry_check, rel=emmo.hasInput)

with SimamsSession() as sess:
    reaxpro_wrapper1 = cuba.Wrapper(session=sess)
    reaxpro_wrapper1.add(calculation,
                        rel=emmo.hasPart)
    reaxpro_wrapper1.session.run()

## Post-processing:
search_calculation = \
        search.find_cuds_objects_by_oclass(
                           emmo.Calculation,
                           reaxpro_wrapper1, emmo.hasPart)

if search_calculation:
    #pretty_print(search_calculation[0])
    Cuds2dot(search_calculation[0]).render()
    export_cuds(search_calculation[0], "CO+Pt111.ttl", format="ttl")

# Binding Site Calculation

calculation = emmo.BindingSites()
num_expeditions = emmo.NumberOfExpeditions(hasNumericalData=1)
num_explorers = emmo.NumberOfExplorers(hasNumericalData=4)

symmetry_check = emmo.CheckSymmetry(hasNumericalData = 'F')
calculation.add(molecule, lattice, solver, fixed_region, num_expeditions,
                forcefield, num_explorers, max_energy, max_distance, ref_region,
                random_seed, symmetry_check, rel=emmo.hasInput)

with SimamsSession() as sess:
    reaxpro_wrapper2 = cuba.Wrapper(session=sess)
    reaxpro_wrapper2.add(calculation,
                        rel=emmo.hasPart)
    reaxpro_wrapper2.session.run()

# Mesoscopic calculation 

calculation = emmo.MesoscopicCalculation()

# pyZacros Mechanism, retrieved from a previous calculation

search_calculation = \
        search.find_cuds_objects_by_oclass(
                           emmo.Calculation,
                           reaxpro_wrapper1, emmo.hasPart)
search_mechanism = \
                   search.find_cuds_objects_by_oclass(
                                           emmo.ChemicalReactionMechanism,
                                           search_calculation[0], emmo.hasOutput)

# pyZacros Cluster_expansions, retrieved from a previous calculation
search_calculation = \
        search.find_cuds_objects_by_oclass(
                           emmo.Calculation,
                           reaxpro_wrapper1, emmo.hasPart)
search_clusters = \
                   search.find_cuds_objects_by_oclass(
                                           emmo.ClusterExpansion,
                                           search_calculation[0], emmo.hasOutput)

for cluster in search_clusters:
    calculation.add(cluster, rel=emmo.hasInput)

# pyZacros Lattice, retrieved from previous calculation, written in file
search_calculation = \
        search.find_cuds_objects_by_oclass(
                           emmo.Calculation,
                           reaxpro_wrapper2, emmo.hasPart)
search_lattice = \
                   search.find_cuds_objects_by_oclass(
                                           crystallography.UnitCell,
                                           search_calculation[0], emmo.hasOutput)

# pyZacros Settings:
#
random_seed = emmo.RandomSeed(hasNumericalData='10')

temperature_float = emmo.Real(hasNumericalData='273.15')
temperature_unit = emmo.Kelvin(hasSymbolData='K')
temperature = emmo.ThermodynamicTemperature()
temperature.add(temperature_unit, rel=emmo.hasReferenceUnit)
temperature.add(temperature_float, rel=emmo.hasQuantityValue)

pressure_float = emmo.Real(hasNumericalData='101325')
pressure_unit = emmo.Pascal(hasSymbolData='Pa')
pressure = emmo.Pressure()
pressure.add(pressure_unit, rel=emmo.hasReferenceUnit)
pressure.add(pressure_float, rel=emmo.hasQuantityValue)

time_float = emmo.Real(hasNumericalData='0.000001')
time_unit = emmo.Second(hasSymbolData='s')
max_time = emmo.MaximumTime()
max_time.add(time_unit, rel=emmo.hasReferenceUnit)
max_time.add(time_float, rel=emmo.hasQuantityValue)

snapshots_float=emmo.Real(hasNumericalData='3.5')
snapshots = emmo.Snapshots(hasSymbolData="on time")
snapshots.add(snapshots_float, rel=emmo.hasSpatialPart)

CO_gas_species = emmo.GasSpecies()
molar_fraction_CO = emmo.AmountFraction()
CO_symbol = emmo.ChemicalElement(hasSymbolData='CO')
molar_fraction_float = emmo.Real(hasNumericalData='0.1')
molar_fraction_CO.add(molar_fraction_float, rel=emmo.hasQuantityValue)
CO_gas_species.add(molar_fraction_CO, CO_symbol, rel=emmo.hasProperty)

calculation.add(search_mechanism[0], search_lattice[0], \
                random_seed, temperature, pressure, max_time, snapshots, \
                CO_gas_species, rel=emmo.hasInput)

with SimzacrosSession() as sess:
    reaxpro_wrapper3 = cuba.Wrapper(session=sess)
    reaxpro_wrapper3.add(calculation,
                        rel=emmo.hasPart)
    reaxpro_wrapper3.session.run()

