from osp.core.namespaces import emmo, cuba, crystallography
from osp.wrappers.simams.simams_session import SimamsSession
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession 
from osp.models.multiscale.co_pt111_meso import COPt111MesoscaleModel
from osp.core.utils import simple_search as search
from osp.core.utils import pretty_print, Cuds2dot, export_cuds
import os

PATH = os.path.dirname(__file__)
molecule = os.path.join(PATH, "XYZ", "CO_ads+Pt111.xyz")
lattice = os.path.join(PATH, "XYZ", "CO_ads+Pt111.xyz")


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
        "random_seed": 10,
        "temperature": 273.15,
        "pressure": 1.01325,
        "n_gas_species": 1,                                                                                                       
        "gas_specs_names": ["CO"],
        "gas_molar_fracs": [0.1],
        "snapshots": ["on time", 3.5],
        "species_numbers": ["on time", 3.5],
        "max_time": 0.00001,
    }
}

model = COPt111MesoscaleModel(**data)

with SimamsSession() as sess:
    reaxpro_wrapper1 = cuba.Wrapper(session=sess)
    reaxpro_wrapper1.add(model.pes_exploration.cuds,
                        rel=emmo.hasPart)
    reaxpro_wrapper1.session.run()

# Binding Site Calculation

with SimamsSession() as sess:
    reaxpro_wrapper2 = cuba.Wrapper(session=sess)
    reaxpro_wrapper2.add(model.binding_site.cuds,
                        rel=emmo.hasPart)
    reaxpro_wrapper2.session.run()

# Mesoscopic calculation 

# pyZacros Mechanism, retrieved from a previous calculation

search_mechanism = reaxpro_wrapper1.get(oclass=emmo.Calculation).pop().get(oclass=emmo.ChemicalReactionMechanism, rel=emmo.hasOutput)

search_clusters = reaxpro_wrapper1.get(oclass=emmo.Calculation).pop().get(oclass=emmo.ClusterExpansion, rel=emmo.hasOutput)

search_lattice = reaxpro_wrapper2.get(oclass=emmo.Calculation).pop().get(oclass=crystallography.UnitCell, rel=emmo.hasOutput)


model.zgb_model.cuds.add(*search_mechanism, *search_lattice, *search_clusters, rel=emmo.hasInput)

with SimzacrosSession() as sess:
    reaxpro_wrapper3 = cuba.Wrapper(session=sess)
    reaxpro_wrapper3.add(model.zgb_model.cuds,
                        rel=emmo.hasPart)
    reaxpro_wrapper3.session.run()

