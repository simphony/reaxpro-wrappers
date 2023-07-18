import scm.pyzacros as pz
from pathlib import Path
from osp.models.zacros.co_pyzacros import COpyZacrosModel 
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession
from osp.core.namespaces import emmo, cuba

EXAMPLES = Path(__file__).absolute().parents[1] / 'tests/test_files'

pz.init()
input_job = pz.ZacrosJob.load_external(EXAMPLES)

pz_job = pz.ZacrosJob(settings=input_job.settings,
                      lattice=input_job.lattice,
                      mechanism=input_job.mechanism,
                      cluster_expansion=input_job.cluster_expansion)

settings_test = {"simulation_input": {
                  "random_seed": 8949321,
                  "temperature": {"value": 900.0},
                  "pressure": 1.0,
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
                }}

lattice_test = {"lattice_input": {"xyz_file": "../tests/test_files/lattice_input.dat"}}


energetics_test = {"energetics_input": [
                   {"name": "CO_Point_brg",
                    "sites": 1,
                    "lattice_state": [1,  "CO*",  1],
                    "site_types": ["brg"],
                    "cluster_eng": -2.36},

                   {"name": "CO_Point_hol",
                    "sites": 1,
                    "lattice_state": [1, "CO*", 1],
                    "site_types": ["hol"],
                    "cluster_eng": -1.85},

                   {"name": "O_Point_brg",
                    "sites": 1,
                    "lattice_state": [1, "O*", 1],
                    "site_types": ["brg"],
                    "cluster_eng": -1.51},

                   {"name": "CO-O_Pair_brg",
                    "sites": 2,
                    "neighboring": ["1-2"],
                    "lattice_state": [[1, "CO*", 1], [2, "O*", 1]],
                    "site_types": ["hol", "hol"],
                    "cluster_eng": 0.05},

                   {"name": "CO2_Bidentat_brg_hol",
                    "sites": 2,
                    "neighboring": ["1-2"],
                    "lattice_state": [[1, "CO2*", 1], [1, "CO2*", 2]],
                    "site_types": ["brg", "hol"],
                    "cluster_eng": -3.42}
                   ]}

mechanism_test = {"mechanism_input": [
                    {"reversible_step": "O2_adsorption",
                     "gas_reacs_prods": ["O2", -1],
                     "sites": 2,
                     "neighboring": ["1-2"],
                     "initial": [[1, "*", 1], [2, "*", 1]],
                     "final": [[1, "O*", 1], [2, "O*", 1]],
                     "variant": {"name": "brg_brg",
                                 "site_types": ["brg", "brg"],
                                 "pre_expon": 7.980e+07,
                                 "pe_ratio":  9.431e-09,
                                 "activ_eng": 0.00}},

                    {"reversible_step": "CO_adsorption",
                     "gas_reacs_prods": ["CO", -1],
                     "sites": 1,
                     "initial": [1, "*", 1],
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
                    ]}

# Syntactic run:
results = pz_job.run()

# Semantic run:
simulation_test = {}
simulation_test.update(**lattice_test, **settings_test, **energetics_test, **mechanism_test)
model = COpyZacrosModel(**simulation_test)
session = SimzacrosSession()
wrapper = cuba.Wrapper(session=session)
wrapper.add(model.cuds, rel=emmo.hasPart)
session.run()


def test_pyZacros_CO_model():

    assert dict(results.get_reaction_network()) == session.get_reaction_network
    assert int(results.number_of_lattice_sites()) == session.number_of_lattice_sites
    assert list(results.gas_species_names()) == session.gas_species_names
    assert list(results.surface_species_names()) == session.surface_species_names
    assert list(results.site_type_names()) == session.site_type_names
    assert int(results.number_of_snapshots()) == session.number_of_snapshots
    assert int(results.number_of_process_statistics()) == session.number_of_process_statistics
    assert list(results.elementary_steps_names()) == session.elementary_steps_names
