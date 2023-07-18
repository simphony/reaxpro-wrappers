"""Run AMS LandscapeRefinment calculation using ReaxPro ontology."""
import os
from osp.core.namespaces import cuba, emmo
from osp.core.utils import Cuds2dot, pretty_print
from osp.core.utils import simple_search as search
from osp.models.ams.energy_landscape_refinement import EnergyLandscapeRefinement
from osp.wrappers.simams.simams_session import SimamsSession

PATH = os.path.dirname(__file__)

XYZ_PATH = os.path.join(PATH, "XYZ")

example = {
    "pathways": [
      {
        "reactant": {
          "xyz_file": os.path.join(XYZ_PATH, "state1.xyz"),
          "charge": -1
        },
        "transition": {
          "xyz_file": os.path.join(XYZ_PATH, "state3.xyz"),
          "charge": -1
        },
        "product": {
          "xyz_file": os.path.join(XYZ_PATH, "state2.xyz"),
          "charge": -1
        }
      },
      {
        "reactant": {
          "xyz_file": os.path.join(XYZ_PATH, "state2.xyz"),
          "charge": -1
        },
        "transition": {
          "xyz_file": os.path.join(XYZ_PATH, "state5.xyz"),
          "charge": -1
        },
        "product": {
          "xyz_file": os.path.join(XYZ_PATH, "state4.xyz"),
          "charge": -1
        }
      }
    ]
  }

model = EnergyLandscapeRefinement(**example)

pretty_print(model.cuds)


with SimamsSession() as sess:
    reaxpro_wrapper = cuba.Wrapper(session=sess)
    reaxpro_wrapper.add(model.cuds, rel=emmo.hasPart)
    reaxpro_wrapper.session.run()
#
## Post-processing:
Cuds2dot(model.cuds).render()
search_calculation = search.find_cuds_objects_by_oclass(
    emmo.Calculation, reaxpro_wrapper, emmo.hasPart
)
if search_calculation:
    pretty_print(search_calculation[0])
