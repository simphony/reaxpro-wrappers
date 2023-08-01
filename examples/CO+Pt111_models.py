from osp.core.namespaces import emmo, cuba, crystallography
from osp.wrappers.simams.simams_session import SimamsSession
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession 
from osp.models.multiscale.co_pt111_meso import COPt111MesoscaleModel
from osp.core.utils import simple_search as search
from osp.core.utils import pretty_print, Cuds2dot, export_cuds
import os

PATH = os.path.dirname(__file__)
molecule = os.path.join("XYZ", "CO_ads+Pt111.xyz")
lattice = os.path.join("XYZ", "CO_ads+Pt111.xyz")
# PES Exploration 

data = {
    "pes_exploration": {
        "molecule": molecule,
        "lattice": lattice,
        "force_field": 'CHONSFPtClNi',
        "solver_type": 'Direct',
        "n_expeditions": 30,
        "n_explorers": 3,
        "max_energy": 2.0,
        "max_distance": 3.8,
        "random_seed": 100, 
        "fixed_region":'surface', 
        "reference_region": 'surface',
        "symmetry_check": 'T', 
        },
    "binding_site": {
        "n_expeditions": 1,
        "n_explorers": 4,
        "symmetry_check": 'F'
    },
    "zgb_model": { 
        "simulation_input": {
            "random_seed": 8949321,
            "temperature": {"value": 900.0},
            "pressure": 1.01325,
            "n_gas_species": 3,                                                                                                       
            "gas_specs_names": ["O2", "CO", "CO2"],
            "gas_energies":  [0.00000e+00, 0.00000e+00, -4.15784e+00],
            "gas_molec_weights": [3.19898e+01, 2.79949e+01, 4.39898e+01],
            "gas_molar_fracs": [5.00000e-01, 5.00000e-01, 0.00000e+00],
            "n_surf_species":   3,  
            "surf_specs_names": ["O*", "CO*", "CO2*"],
            "surf_specs_dent": [1, 1, 2],
            "snapshots": ["on time", 0.0005],
            "process_statistics": ["on time", 0.0005],
            "species_numbers": ["on time", 0.0005],
            "event_report": "off",
            "max_steps": "infinity",
            "max_time": 200.0,
            "wall_time": 100 
        },

    "lattice_input": {"xyz_file": "XYZ/lattice_input.dat"},

    "energetics_input": [
        { "name": "CO_Point_brg",
        "sites": 1,
        "lattice_state": [1,  "CO*",  1], 
        "site_types": ["brg"],
        "cluster_eng": -2.36},
            
        { "name": "CO_Point_hol",
        "sites": 1,
        "lattice_state": [1,  "CO*",  1], 
        "site_types": ["hol"],
        "cluster_eng": -1.85},
            
        { "name": "O_Point_brg",
        "sites": 1,
        "lattice_state": [1,  "O*",  1], 
        "site_types": ["brg"],
        "cluster_eng": -1.51},
            
        { "name": "CO-O_Pair_brg",
        "sites": 2,
        "neighboring": ["1-2"],
        "lattice_state": [[1,  "CO*",  1], [2, "O*", 1] ],
        "site_types": ["hol", "hol"],
        "cluster_eng": 0.05},
            
        { "name": "CO2_Bidentat_brg_hol",
        "sites": 2,
        "neighboring": ["1-2"],
        "lattice_state": [[1,  "CO2*",  1], [1, "CO2*", 2] ],
        "site_types": ["brg", "hol"],
        "cluster_eng": -3.42}
        ], 

    "mechanism_input": [
            {"reversible_step": "O2_adsorption",
            "gas_reacs_prods": ["O2", -1],
            "sites": 2,
            "neighboring": ["1-2"],
            "initial": [[ 1, "*", 1], [ 2, "*", 1]],
            "final": [[1, "O*", 1], [2, "O*", 1]],
            "variant": {"name": "brg_brg",
                        "site_types": ["brg", "brg"],
                        "pre_expon": 7.980e+07,
                        "pe_ratio":  9.431e-09,
                        "activ_eng": 0.00}},

            {"reversible_step": "CO_adsorption",
            "gas_reacs_prods": ["CO", -1],
            "sites": 1,
            "initial": [ 1, "*", 1],
            "final":   [1, "CO*", 1],
            "variant": {"name": "brg",
                        "site_types": ["brg"],
                        "pre_expon": 4.265e+07,
                        "pe_ratio":  6.563e-09,
                        "activ_eng": 0.00}},

            {"reversible_step": "CO_O_oxidation",
            "gas_reacs_prods": ["CO2", 1],
            "sites": 2,
            "neighboring": ["1-2"],
            "initial": [[1, "CO*", 1], [2, "O*", 1]],
            "final":   [[1, "*", 1], [1, "*", 1]],
            "variant": {"name": "brg_brg",
                        "site_types": ["brg", "brg"],
                        "pre_expon": 2.786e+12,
                        "pe_ratio":  3.231e+07,
                        "activ_eng": 0.52}}
    ]
    }

}

model = COPt111MesoscaleModel(**data)

with SimamsSession() as sess:
    reaxpro_wrapper1 = cuba.Wrapper(session=sess)
    reaxpro_wrapper1.add(model.pes_exploration.cuds,
                        rel=emmo.hasPart)
    reaxpro_wrapper1.session.run()

## Post-processing:
search_calculation = \
        search.find_cuds_objects_by_oclass(
                           emmo.Calculation,
                           reaxpro_wrapper1, emmo.hasPart)

if search_calculation:
    #pretty_print(search_calculation[0])
    Cuds2dot(model.pes_exploration.cuds).render()
    export_cuds(model.pes_exploration.cuds, "CO+Pt111.ttl", format="ttl")

# Binding Site Calculation

with SimamsSession() as sess:
    reaxpro_wrapper2 = cuba.Wrapper(session=sess)
    reaxpro_wrapper2.add(model.binding_site.cuds,
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
snapshots = emmo.Snapshots(hasSymbolData="off")
snapshots.add(snapshots_float, rel=emmo.hasSpatialPart)

molar_fraction_CO = emmo.AmountFraction()
CO_symbol = emmo.ChemicalElement(hasSymbolData='CO')
molar_fraction_CO.add(CO_symbol, rel=emmo.hasProperty)
molar_fraction_float = emmo.Real(hasNumericalData='0.1')
molar_fraction_CO.add(molar_fraction_float, rel=emmo.hasQuantityValue)

calculation.add(search_mechanism[0], search_lattice[0], \
                random_seed, temperature, pressure, max_time, snapshots, \
                molar_fraction_CO, rel=emmo.hasInput)

with SimzacrosSession() as sess:
    reaxpro_wrapper3 = cuba.Wrapper(session=sess)
    reaxpro_wrapper3.add(calculation,
                        rel=emmo.hasPart)
    reaxpro_wrapper3.session.run()

