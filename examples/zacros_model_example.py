"""Example of py/Zacros model"""

import os
import json
from osp.models.zacros.co_pyzacros import COpyZacrosModel
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession
from osp.core.namespaces import emmo, cuba
# from osp.core.utils import pretty_print

#path = os.path.join(os.path.dirname(__file__), "Ziff-Gulari-Barshad-model.json")
path = os.path.join(os.path.dirname(__file__), "Ziff-Gulari-Barshad-model-variants.json")

with open(path, mode="r+") as file:
    content = json.loads(file.read())


model = COpyZacrosModel(**content)
# pretty_print(model.cuds)


session = SimzacrosSession()
wrapper = cuba.Wrapper(session=session)
wrapper.add(model.cuds, rel=emmo.hasPart)
session.run()
# print(".json content:", content)
