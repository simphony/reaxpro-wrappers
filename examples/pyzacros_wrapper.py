"""ReaxPro example, adding a semantic mechanism and cluster_expansions to a PyZacros/Zacros calculation."""

from osp.core.namespaces import emmo, cuba
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession 
from osp.tools.io_functions import read_mechanism, read_cluster_expansion
from osp.core.utils import pretty_print
import scm.plams

scm.plams.JobRunner(parallel=True, maxjobs=4)

calculation = emmo.MesoscopicCalculation()

# Add Mechanism from file:
mechanism=read_mechanism('./Zacros_inputs/mechanism_input.dat')
calculation.add(mechanism, rel=emmo.hasInput)

# Add ClusterExpansions from file:
clusters=read_cluster_expansion('./Zacros_inputs/energetics_input.dat')

for cluster in clusters:
    calculation.add(cluster, rel=emmo.hasInput)

# Add different Settings of the calculation:
random_seed = emmo.RandomSeed(hasNumericalData='8949321')

temperature_float = emmo.Real(hasNumericalData='900.0')
temperature_unit = emmo.Kelvin(hasSymbolData='K')
temperature = emmo.ThermodynamicTemperature()
temperature.add(temperature_unit, rel=emmo.hasReferenceUnit)
temperature.add(temperature_float, rel=emmo.hasQuantityValue)

pressure_float = emmo.Real(hasNumericalData='100000')
pressure_unit = emmo.Pascal(hasSymbolData='Pa')
pressure = emmo.Pressure()
pressure.add(pressure_unit, rel=emmo.hasReferenceUnit)
pressure.add(pressure_float, rel=emmo.hasQuantityValue)

time_float = emmo.Real(hasNumericalData='200')
time_unit = emmo.Second(hasSymbolData='s')
max_time = emmo.MaximumTime()
max_time.add(time_unit, rel=emmo.hasReferenceUnit)
max_time.add(time_float, rel=emmo.hasQuantityValue)

time_float = emmo.Real(hasNumericalData='100')
time_unit = emmo.Second(hasSymbolData='s')
wall_time = emmo.WallTime()
wall_time.add(time_unit, rel=emmo.hasReferenceUnit)
wall_time.add(time_float, rel=emmo.hasQuantityValue)

max_steps_float=emmo.Real(hasNumericalData='9223372036854775807')
max_steps = emmo.MaximumSteps()
max_steps.add(max_steps_float, rel=emmo.hasSpatialPart)

snapshots_float=emmo.Real(hasNumericalData='0.0005')
snapshots = emmo.Snapshots(hasSymbolData="on time")
snapshots.add(snapshots_float, rel=emmo.hasSpatialPart)

molar_fraction_CO = emmo.AmountFraction()
CO_symbol = emmo.ChemicalElement(hasSymbolData='CO')
molar_fraction_CO.add(CO_symbol, rel=emmo.hasProperty)
molar_fraction_float = emmo.Real(hasNumericalData='0.5')
molar_fraction_CO.add(molar_fraction_float, rel=emmo.hasQuantityValue)

molar_fraction_O2 = emmo.AmountFraction()
O2_symbol = emmo.ChemicalElement(hasSymbolData='O2')
molar_fraction_O2.add(O2_symbol, rel=emmo.hasProperty)
molar_fraction_float = emmo.Real(hasNumericalData='0.5')
molar_fraction_O2.add(molar_fraction_float, rel=emmo.hasQuantityValue)

calculation.add(random_seed, temperature, pressure, \
                max_time, wall_time, max_steps,\
                snapshots, \
                molar_fraction_CO, molar_fraction_O2, rel=emmo.hasInput)

#pretty_print(calculation)

# Run the calculation:
with SimzacrosSession() as sess:
    reaxpro_wrapper = cuba.Wrapper(session=sess)
    reaxpro_wrapper.add(calculation,
                        rel=emmo.hasPart)
    reaxpro_wrapper.session.run()
