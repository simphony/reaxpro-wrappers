# AMS run using the default XC_functional, basis_set and calculation_type:

from osp.core.namespaces import emmo, cuba
from osp.tools.io_functions import read_molecule
from osp.wrappers.simams.simams_session import SimamsSession
from osp.core.utils import pretty_print
from osp.core.utils import Cuds2dot
from osp.core.utils import import_cuds, export_cuds

molecule = read_molecule('./XYZ/H2O.xyz')
calculation = emmo.Calculation()
calculation.add(molecule, rel=emmo.hasInput)

with SimamsSession() as sess:
    reaxpro_wrapper = cuba.Wrapper(session=sess)
    reaxpro_wrapper.add(calculation, rel=emmo.hasPart)
    reaxpro_wrapper.session.run()
# Post-processing:
pretty_print(reaxpro_wrapper)
Cuds2dot(reaxpro_wrapper).render()
export_cuds(reaxpro_wrapper,"ams_wrapper.ttl",format="ttl")           
