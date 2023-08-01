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
        "random_seed": 10,
        "temperature": 273.15,
        "pressure": 101325,
        "n_gas_species": 1,                                                                                                       
        "gas_specs_names": ["CO"],
        "gas_molar_fracs": [0.1],
        "snapshots": ["on time", 3.5],
        "max_time": 0.000001,
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


# pyZacros Mechanism, retrieved from a previous calculation

search_mechanism = \
                   search.find_cuds_objects_by_oclass(
                                           emmo.ChemicalReactionMechanism,
                                           model.pes_exploration.cuds, emmo.hasOutput)

# pyZacros Cluster_expansions, retrieved from a previous calculation

search_clusters = \
                   search.find_cuds_objects_by_oclass(
                                           emmo.ClusterExpansion,
                                           model.pes_exploration.cuds, emmo.hasOutput)

for cluster in search_clusters:
    model.zgb_model.cuds.add(cluster, rel=emmo.hasInput)

# pyZacros Lattice, retrieved from previous calculation, written in file

search_lattice = \
                   search.find_cuds_objects_by_oclass(
                                           crystallography.UnitCell,
                                           model.binding_site.cuds, emmo.hasOutput)

# pyZacros Settings:
#

model.zgb_model.cuds.add(search_mechanism[0], search_lattice[0], rel=emmo.hasInput)

with SimzacrosSession() as sess:
    reaxpro_wrapper3 = cuba.Wrapper(session=sess)
    reaxpro_wrapper3.add(model.zgb_model.cuds,
                        rel=emmo.hasPart)
    reaxpro_wrapper3.session.run()

