import os
from osp.core.namespaces import emmo, cuba, crystallography
from osp.wrappers.simams.simams_session import SimamsSession
from osp.wrappers.simzacros.simzacros_session import SimzacrosSession 
from osp.models.multiscale.co_pt111_meso import COPt111MesoscaleModel

# Mesoscopic calculation
from osp.core.utils import import_cuds


with SimzacrosSession() as sess:
    reaxpro_wrapper2 = cuba.Wrapper(session=sess)
    imports = import_cuds("export.ttl", session = sess)
    if isinstance(imports, list):

        reaxpro_wrapper2.add(*imports,
                             rel=cuba.relationship)
    else:
        reaxpro_wrapper2.add(imports,
                             rel=cuba.relationship)
    reaxpro_wrapper2.session.run()

from osp.core.utils import pretty_print

workflow = reaxpro_wrapper2.get(rel=cuba.relationship).pop()
pretty_print(workflow)
