import os
from osp.core.namespaces import emmo, cuba, crystallography
from osp.wrappers.simams.simams_session import SimamsSession
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession 
from osp.models.multiscale.co_pt111_meso import COPt111MesoscaleModel

PATH = os.path.dirname(__file__)
molecule = os.path.join(PATH, "XYZ", "CO_ads+Pt111.xyz")
lattice = os.path.join(PATH, "XYZ", "CO_ads+Pt111.xyz")

# setup data model

data = {
    "pes_exploration": {
        "molecule": molecule,
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
        "n_expeditions": 30,
        "n_explorers": 4,
        "symmetry_check": 'F'
    },
    "zgb_model": {
        "random_seed": 10,
        "temperature": 273.15,
        "pressure": 1.01325,
        "n_gas_species": 1,                                                                                                       
        "gas_specs_names": ["CO"],
        "gas_molar_fracs": [0.1],
        "snapshots": ["on time", 3.5],
        "species_numbers": ["on time", 3.5],
        "process_statistics": ["on time", 3.5],
        "max_time": 0.00001,
    },
    "adp": {
        "min": 0.2,
        "max": 0.8,
        "num": 5
    }
}

model = COPt111MesoscaleModel(**data)

# PES Exploration + Binding Site calculation

with SimamsSession() as sess:
    reaxpro_wrapper1 = cuba.Wrapper(session=sess)
    reaxpro_wrapper1.add(model.cuds,
                        rel=cuba.relationship)
    reaxpro_wrapper1.session.run()


# map output from previous calculation to next calculation

workflow = reaxpro_wrapper1.get(rel=cuba.relationship).pop()

process_search = workflow.get(oclass=emmo.ProcessSearch).pop()

binding_sites = workflow.get(oclass=emmo.BindingSites).pop()

mesocopic = workflow.get(oclass=emmo.MesoscopicCalculation).pop()

search_mechanism = process_search.get(oclass=emmo.ChemicalReactionMechanism, rel=emmo.hasOutput)

search_clusters =  process_search.get(oclass=emmo.ClusterExpansion, rel=emmo.hasOutput)

search_lattice = binding_sites.get(oclass=crystallography.UnitCell, rel=emmo.hasOutput)

mesocopic.add(*search_mechanism, *search_lattice, *search_clusters, rel=emmo.hasInput)

# Mesoscopic calculation

with SimzacrosSession() as sess:
    reaxpro_wrapper2 = cuba.Wrapper(session=sess)
    reaxpro_wrapper2.add(workflow,
                         rel=cuba.relationship)
    reaxpro_wrapper2.session.run()

from osp.core.utils import pretty_print

workflow = reaxpro_wrapper2.get(rel=cuba.relationship).pop()
pretty_print(workflow)
