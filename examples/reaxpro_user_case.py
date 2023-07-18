"""ReaxPro platform test case.

This sematic-based script showcase the connection with MarketPlaces.
It connects two engines AMS and Zacros via their syntactic engines
(PLAMS and pyZacros) using different interoperability layers.

This script has no scientific meaning per se. The mission here
is to make a dummy connection workflow. The main steps are:

1. Calculate CO2 formation energy via PLAMS/AMS.
2. Upload the float value to a MarketPlace.
3. Retrieve the value from the MarketPlace.
4. Print out simulation_input.dat with the energy value and perform
   the KMC simulation.
"""

from pathlib import Path

import pyzacros as pz
import scm.plams
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession
#from scm.plams import init as PlamsInit
#
from osp.core.namespaces import cuba, emmo
from osp.wrappers.simams.simams_session import SimamsSession
from osp.core.utils import pretty_print
from osp.tools.io_functions import read_molecule
from osp.tools.mapping_functions import map_energies
# from osp.core.utils import serialize_rdf_graph,import_rdf_file

# Uncomment here to use a energy value and avoid the long-running block below
#Formation_energy_CO2_value = -4.157842190497977
#
#  1. Calculate CO2 formation energy via PLAMS/AMS.
simulation_AMS = emmo.Simulation()
CO_gas = read_molecule('./XYZ/CO.xyz')
CO2_gas = read_molecule('./XYZ/CO2.xyz')
O2_gas = read_molecule('./XYZ/O2.xyz')
basis_set = emmo.DZP()
xc_functional = emmo.PBE()

molecules = [CO_gas, CO2_gas, O2_gas]  # all of them are emmo.Molecule() objects
geom_opt_calculations = [emmo.GeometryOptimization() for i in range(3)]

for index, molecule in enumerate(molecules):
    geom_opt_calculations[index].add(molecule,
                                     basis_set,
                                     xc_functional,
                                     rel=emmo.hasPart)
simulation_AMS.add(geom_opt_calculations[0], rel=emmo.hasTemporalFirst)
simulation_AMS.add(geom_opt_calculations[1], rel=emmo.hasTemporalNext)
simulation_AMS.add(geom_opt_calculations[2], rel=emmo.hasTemporalLast)

#   Run AMS
with SimamsSession() as sess:
    reaxpro_wrapper = cuba.Wrapper(session=sess)
    reaxpro_wrapper.add(simulation_AMS, rel=emmo.hasPart)
    reaxpro_wrapper.session.run()

##    Retrieve total energies:
list_of_energies = map_energies(reaxpro_wrapper)
E_CO = float(list_of_energies[0])
E_CO2 = float(list_of_energies[1])
E_O2 = float(list_of_energies[2])
#    Formation_energy_CO2 = E_CO2 - (E_CO + 1/2*E_O2)
Formation_energy_CO2_value = E_CO2 - (E_CO + 1 / 2 * E_O2)
print('Formation energy CO2:', Formation_energy_CO2_value, 'eV')
Formation_energy_CO2 = emmo.FormationEnergy()
energy_value = emmo.Real(hasNumericalData=Formation_energy_CO2_value)
energy_unit = emmo.ElectronVolt(hasSymbolData='eV')
Formation_energy_CO2.add(energy_value, rel=emmo.hasQuantityValue)
Formation_energy_CO2.add(energy_unit, rel=emmo.hasReferenceUnit)

#   2. Upload the float value to a MarketPlace.

# export_cuds(Formation_energy_CO2, "Formation_energy_CO2.ttl", format="ttl")
# print("pretty_print before the database upload:")
# HOST = '127.0.0.1'  # It needs to be filled with the host IP
# PORT = 8687         # and port.
# DB = 'reaxpro_calculations.db'
# with TransportSessionClient(SqliteSession, uri='ws://188.166.162.208:8687',
#                            path=DB) as session:
#    wrapper = cuba.Wrapper(session=session)
#    wrapper.add(Formation_energy_CO2, rel=emmo.hasPart)
#    # Add the CUDs object to the database
#    session.commit()
#
#   3. Retrieve the value from the MarketPlace.
# with TransportSessionClient(SqliteSession, uri='ws://188.166.162.208:8687',
#                            path=DB) as session:
#    wrapper = cuba.Wrapper(session=session)
#    energy = wrapper.get(oclass=emmo.TotalElectronicEnergy)[0]
#    print("pretty_print after the database upload/dowload:")
#    pretty_print(energy)
#     wrapper = cuba.wrapper(session=session)
#
#   4. Print out simulation_input.dat with the energy value and perform
#      the KMC simulation.

#
# Semantics
s0 = emmo.AdsorbedSpecies()
s0_integer = emmo.Integer(hasNumericalData=1)
s0_denticity = emmo.Denticity()
s0_denticity.add(s0_integer, rel=emmo.hasQuantityValue)
s0.add(s0_denticity, rel=emmo.hasProperty)

O2_gas = emmo.GasSpecies()
O2_gas_symbol = emmo.ChemicalElement(hasSymbolData='O2')

O_adsorbed = emmo.AdsorbedSpecies()
O_adsorbed_integer = emmo.Integer(hasNumericalData=1)
O_adsorbed_denticity = emmo.Denticity()
O_adsorbed_denticity.add(O_adsorbed_integer, rel=emmo.hasQuantityValue)
O_adsorbed_symbol = emmo.ChemicalElement(hasSymbolData='O')
O_adsorbed.add(O_adsorbed_symbol, O_adsorbed_denticity, rel=emmo.hasProperty)

CO_adsorbed = emmo.AdsorbedSpecies()
CO_adsorbed_integer = emmo.Integer(hasNumericalData=1)
CO_adsorbed_denticity = emmo.Denticity()
CO_adsorbed_denticity.add(CO_adsorbed_integer, rel=emmo.hasQuantityValue)
CO_adsorbed_symbol = emmo.ChemicalElement(hasSymbolData='CO')
CO_adsorbed.add(CO_adsorbed_symbol,
                CO_adsorbed_denticity,
                rel=emmo.hasProperty)

CO2_adsorbed = emmo.AdsorbedSpecies()
CO2_adsorbed_integer = emmo.Integer(hasNumericalData=2)
CO2_adsorbed_denticity = emmo.Denticity()
CO2_adsorbed_denticity.add(CO2_adsorbed_integer, rel=emmo.hasQuantityValue)
CO2_adsorbed_symbol = emmo.ChemicalElement(hasSymbolData='CO2')
CO2_adsorbed.add(CO2_adsorbed_symbol,
                 CO2_adsorbed_denticity,
                 rel=emmo.hasProperty)

calculation_KMC = emmo.calculation()

random_seed = emmo.SimulationSetting()
random_seed_integer = emmo.Integer(hasNumericalData=8949321)
random_seed.add(random_seed_integer, rel=emmo.hasQuantityValue)
pretty_print(random_seed)

KMC_engine = emmo.Zacros()

temperature_float = emmo.Real(hasNumericalData=900.0)
temperature_unit = emmo.Kelvin(hasSymbolData='K')
temperature = emmo.ThermodynamicTemperature()
temperature.add(temperature_float, rel=emmo.hasQuantityValue)
pretty_print(temperature)
# export_cuds(test_boolean,"reaxpro_test.ttl",format="ttl")

pressure_float = emmo.Real(hasNumericalData=100000)
pressure_unit = emmo.Pascal(hasSymbolData='Pa')
pressure = emmo.Pressure()
pressure.add(pressure_unit, rel=emmo.hasReferenceUnit)
pressure.add(pressure_float, rel=emmo.hasQuantityValue)
pretty_print(pressure)

wall_time_float = emmo.Real(hasNumericalData=30)
wall_time_unit = emmo.Second(hasSymbolData='s')
wall_time = emmo.WallTime()
wall_time.add(wall_time_unit, rel=emmo.hasReferenceUnit)
wall_time.add(wall_time_float, rel=emmo.hasQuantityValue)
#pretty_print(wall_time)
#
## End Semantics
#
## Instantiate settings:
sett = pz.Settings()
# Species:
s0 = pz.Species('*', 1)  # Empty adsorption site
# - Gas-species:
O2_gas = pz.Species('O2', gas_energy=0.0)
sett.molar_fraction.O2 = 0.50
CO_gas = pz.Species('CO')
sett.molar_fraction.CO = 0.50
#CO2_gas = Species('CO2', gas_energy=Formation_energy_CO2_value)
CO2_gas = pz.Species("CO2", gas_energy=-4.15784) # original value
#
# Adsorbed species:
O_adsorbed = pz.Species('O*', 1)
CO_adsorbed = pz.Species('CO*', 1)
CO2_adsorbed = pz.Species('CO2*', 2)
#
# O2_adsorption:
O2_adsorption = pz.ElementaryReaction(site_types=['brg', 'brg'],
                                   initial=[s0, s0, O2_gas],
                                   final=[O_adsorbed, O_adsorbed],
                                   neighboring=[(0, 1)],
                                   reversible=True,
                                   pre_expon=7.980e+07,
                                   pe_ratio=9.431e-09,
                                   activation_energy=0.0)

# CO_adsoprtion:
CO_adsorption = pz.ElementaryReaction(site_types=['brg'],
                                   initial=[s0, CO_gas],
                                   final=[CO_adsorbed],
                                   reversible=True,
                                   pre_expon=4.265e+07,
                                   pe_ratio=6.563e-09,
                                   activation_energy=0.0)

# CO_oxidation:
CO_oxidation = pz.ElementaryReaction(site_types=['brg', 'brg'],
                                  initial=[CO_adsorbed, O_adsorbed],
                                  final=[s0, s0, CO2_gas],
                                  neighboring=[(0, 1)],
                                  reversible=True,
                                  pre_expon=2.786e+12,
                                  pe_ratio=3.231e+07,
                                  activation_energy=0.52)

# Settings:
sett.random_seed = 8949321
sett.temperature = 900.0
sett.pressure = 1.00
sett.snapshots = ('time', 5.e-4)
sett.process_statistics = ('time', 5.e-4)
sett.species_numbers = ('time', 5.e-4)
sett.event_report = 'off'
sett.max_steps = 'infinity'
sett.max_time = 200.0
sett.wall_time = 100
#
inputs = Path(__file__).parent / 'Zacros_inputs'
input_job = pz.ZacrosJob.load_external(inputs)
## Syntactic 
myJob = pz.ZacrosJob(
    settings=sett,
#    mechanism=[O2_adsorption, CO_adsorption, CO_oxidation],
    mechanism=input_job.mechanism,
    lattice=input_job.lattice,
    cluster_expansion=input_job.cluster_expansion)

pz.init()
results = myJob.run()

## Run Zacros using semantics:
#calculation_KMC = emmo.MesoscopicCalculation()
#calculation_KMC2 = emmo.MesoscopicCalculation()
#with SimzacrosSession() as sess:
#    reaxpro_wrapper = cuba.Wrapper(session=sess)
#    reaxpro_wrapper.add(calculation_KMC, calculation_KMC2,rel=emmo.hasPart)
#    reaxpro_wrapper.session.run()
#
