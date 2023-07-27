"""Example of py/Zacros model"""

import os
import json
from osp.models.zacros.co_pyzacros import COpyZacrosModel
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession
from osp.core.namespaces import emmo, cuba
from osp.core.utils import import_cuds
from osp.core.utils import pretty_print

#path = os.path.join(os.path.dirname(__file__), "Ziff-Gulari-Barshad-model.json")
path = os.path.join(os.path.dirname(__file__), "Ziff-Gulari-Barshad-model-variants.json")

with open(path, mode="r+") as file:
    content = json.loads(file.read())

# alternatively with the standard file. Needs connection to minio!
#content = COpyZacrosModel.Config.schema_extra["example"]
model = COpyZacrosModel(**content)


session = SimzacrosSession()
wrapper = cuba.Wrapper(session=session)
cuds = import_cuds(model.file, session=session)
wrapper.add(*cuds, rel=cuba.relationship)
session.run()
# print(".json content:", content)
