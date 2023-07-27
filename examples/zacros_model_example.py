"""Example of py/Zacros model"""

import os
from osp.models.zacros.co_pyzacros import COpyZacrosModel
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession
from osp.core.namespaces import cuba
from osp.core.utils import import_cuds


# alternatively with the standard file. Needs connection to minio!
content = COpyZacrosModel.Config.schema_extra["example"]
directory = os.path.dirname(__file__)
content["lattice_input"]["xyz_file"] = os.path.join(directory, "XYZ", "lattice_input.dat")
model = COpyZacrosModel(**content)


session = SimzacrosSession()
wrapper = cuba.Wrapper(session=session)
cuds = import_cuds(model.file, session=session)
wrapper.add(*cuds, rel=cuba.relationship)
session.run()
# print(".json content:", content)
