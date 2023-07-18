"""Run ADF/PLAMS using the ReaxPro ontology."""

from osp.core.namespaces import emmo, cuba
from osp.tools.io_functions import read_molecule
from osp.wrappers.simams.simams_session import SimamsSession
from osp.core.utils import simple_search as search
from osp.core.utils import pretty_print 
from os import listdir


path_with_geometries = "XYZ"
energies = []
molecules = listdir(path_with_geometries)
for i in molecules:
    molecule = read_molecule(path_with_geometries+"/"+i)
    basis_set = emmo.DZP()
    xc_functional = emmo.PBE()
    calculation = emmo.GeometryOptimization()
    calculation.add(molecule, basis_set, xc_functional,
                    rel=emmo.hasPart)

    with SimamsSession() as sess:
        reaxpro_wrapper = cuba.Wrapper(session=sess)
        reaxpro_wrapper.add(calculation, rel=emmo.hasPart)
        reaxpro_wrapper.session.run()
        list_empty = search.find_cuds_objects_by_oclass(
                           oclass=emmo.Real,
                           root=reaxpro_wrapper,
                           rel=emmo.hasPart)
        energies.append(list_empty[0].hasNumericalData)

print("Total electronic energies of optimized structures, in eV:")
print(list(zip(molecules, energies)))
