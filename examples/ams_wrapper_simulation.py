"""Run a dummy AMS calculation using a Simulation CUDS"""

from tkinter import W
from osp.core.namespaces import emmo, cuba
from osp.tools.io_functions import read_molecule
from osp.tools.mapping_functions import map_energies 
from osp.wrappers.simams.simams_session import SimamsSession
from osp.core.utils import pretty_print
from osp.core.utils import simple_search as search

simulation = emmo.Simulation()

molecule1 = read_molecule('./XYZ/H2O.xyz')
basis_set1 = emmo.SZ()
xc_functional1 = emmo.LDA()
calculation1 = emmo.GeometryOptimization()

molecule2 = read_molecule('./XYZ/H2O.xyz')
basis_set2 = emmo.DZP()
xc_functional2 = emmo.PBE()
calculation2 = emmo.WavefunctionOptimization()

molecule3 = read_molecule('./XYZ/H2O.xyz')
basis_set3 = emmo.TZP()
xc_functional3 = emmo.B3LYP()
calculation3 = emmo.GeometryOptimization()

calculation1.add(molecule1, basis_set1,
                 xc_functional1,
                 rel=emmo.hasInput)

calculation2.add(molecule2, basis_set2,
                 xc_functional2,
                 rel=emmo.hasInput)

calculation3.add(molecule3, basis_set3,
                 xc_functional3,
                 rel=emmo.hasInput)

simulation.add(calculation1, rel=emmo.hasTemporalFirst)
simulation.add(calculation2, rel=emmo.hasTemporalNext)
simulation.add(calculation3, rel=emmo.hasTemporalLast)

with SimamsSession() as sess:
    reaxpro_wrapper = cuba.Wrapper(session=sess)
    reaxpro_wrapper.add(simulation,
                        rel=emmo.hasPart)
    reaxpro_wrapper.session.run()

print ("List of total energies:", map_energies(reaxpro_wrapper), "eV")
#pretty_print(reaxpro_wrapper)