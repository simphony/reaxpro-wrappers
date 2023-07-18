"""Run AMS calculation using ReaxPro ontology.""" 

from osp.core.namespaces import emmo, cuba
from osp.tools.io_functions import read_molecule
from osp.wrappers.simams.simams_session import SimamsSession
from osp.core.utils import pretty_print
from osp.core.utils import Cuds2dot
from osp.core.utils import export_cuds

# Read molecule from xyz format:
#molecule = read_molecule('./XYZ/H2O.xyz')
# Read molecule from mol format:
molecule = read_molecule('./XYZ/ChEBI_15377.mol')

pretty_print(molecule)
basis_set = emmo.DZP()
xc_functional = emmo.PBE()
calculation = emmo.GeometryOptimization()

calculation.add(molecule, basis_set,
                xc_functional,
                rel=emmo.hasInput)

with SimamsSession() as sess:
    reaxpro_wrapper = cuba.Wrapper(session=sess)
    reaxpro_wrapper.add(calculation,
                        rel=emmo.hasPart)
    reaxpro_wrapper.session.run()

# Post-processing:
pretty_print(reaxpro_wrapper)
Cuds2dot(reaxpro_wrapper).render()
export_cuds(reaxpro_wrapper, "ams_wrapper.ttl", format="ttl")
