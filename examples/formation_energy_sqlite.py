"""Run ADF/PLAMS using the co2_activation ontology."""

from osp.core.namespaces import emmo, cuba
from osp.tools.io_functions import read_molecule
from osp.wrappers.simams.simams_session import SimamsSession
from osp.core.utils import simple_search as search
from os import listdir
from osp.core.utils import pretty_print
from osp.wrappers.sqlite.sqlite_session import SqliteSession
from osp.core.session.transport.transport_session_client import \
    TransportSessionClient


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
        reaxpro_wrapper.add(molecule, basis_set,
                            xc_functional, calculation,
                            rel=emmo.hasPart)
        reaxpro_wrapper.session.run()
        list_empty = search.find_cuds_objects_by_oclass(
                           oclass=emmo.Real,
                           root=reaxpro_wrapper,
                           rel=emmo.hasPart)
        energy = search.find_cuds_objects_by_oclass(
                           oclass=emmo.TotalElectronicEnergy,
                           root=reaxpro_wrapper,
                           rel=emmo.hasPart)
        energies.append(list_empty[0].hasNumericalData)

print("Total electronic energies of optimized structures, in eV:")
print(list(zip(molecules, energies)))

print("pretty_print before the database upload:")
pretty_print(energy[0])
HOST = '127.0.0.1'  # It needs to be filled with the host IP
PORT = 8687         # and port.
DB = 'reaxpro_calculations.db'
with TransportSessionClient(SqliteSession, uri='ws://188.166.162.208:8687',
                            path=DB) as session:
    wrapper = cuba.Wrapper(session=session)
    wrapper.add(energy[0], rel=emmo.hasPart)
    # Add the CUDs object to the database
    session.commit()

with TransportSessionClient(SqliteSession, uri='ws://188.166.162.208:8687',
                            path=DB) as session:
    wrapper = cuba.Wrapper(session=session)
    energy = wrapper.get(oclass=emmo.TotalElectronicEnergy)[0]
    print("pretty_print after the database upload/dowload:")
    pretty_print(energy)
#     wrapper = cuba.wrapper(session=session)
