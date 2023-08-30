"""Example of py/Zacros model"""

import os
from osp.models.zacros.co_pyzacros import COpyZacrosModel
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession
from osp.core.namespaces import cuba


# alternatively with the standard file. Needs connection to minio!
content = COpyZacrosModel.Config.schema_extra["example"]
directory = os.path.dirname(__file__)
content["lattice_input"]["file"] = os.path.join(directory, "XYZ", "lattice_input.dat")
model = COpyZacrosModel(**content)

with SimzacrosSession() as session:
    wrapper = cuba.Wrapper(session=session)
    wrapper.add(model.cuds, rel=cuba.relationship)
    session.run()
# print(".json content:", content)
