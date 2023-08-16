"""Describe AMS Session class."""
from osp.core.session import SimWrapperSession
from osp.core.cuds import Cuds
from osp.tools.graph_functions import graph_wrapper_dependencies
from osp.tools.mapping_functions import map_function, map_results
from scm.plams import MultiJob, AMSJob, init, config, JobRunner
from uuid import uuid4
import os
import multiprocessing
# from osp.core.utils import pretty_print


class SimamsSession(SimWrapperSession):
    """Describe AMS Session class."""

    def __init__(self, engine=None, **kwargs):
        """Initialise SimamsSession."""
        if engine is None:
            init()
            maxjobs = multiprocessing.cpu_count()
            n_procs = os.environ.get("REAXPRO_N_PROCESSES") or 1
            if isinstance(n_procs, str):
                n_procs = int(n_procs)
            config.default_jobrunner = JobRunner(parallel=True, maxjobs=maxjobs)
            config.job.runscript.nproc = n_procs
            self.workdir = config.get("default_jobmanager").workdir
            self.jobname = str(uuid4())
            self.engine = MultiJob(name=self.jobname)
        super().__init__(engine)

    def __str__(self):
        """To overwrite the private str method. Not advised, but here it is."""
        # TODO: Define the output of str(SomeSimulationSession())
        return "Some Wrapper Session"

    def _run(self, root_cuds_object: Cuds):
        """
        Execute engine session.

        :param root_cuds_object: Main Cuds object.

        """
        self.engine.run()
        self._tarball = map_results(self.engine, root_cuds_object)

    # OVERRIDE   # Map results
    def _load_from_backend(self, uids, expired=None):
        """Load the cuds object from the simulation engine."""
        # TODO load cuds objects from the backend
        for uid in uids:
            if uid in self._registry:
                yield self._registry.get(uid)
            else:
                yield None

    # OVERRIDE #Map functions
    def _apply_added(self, root_obj, buffer):
        """Add the added cuds to the engine."""
        # Solve dependencies:
        self.dependencies = graph_wrapper_dependencies(root_obj)

        # One Calculation without Simulation:
        if list(self.dependencies.keys())[0] == "Calculation":
            calculation = self.dependencies["Calculation"][0]
            (plams_molecule, plams_settings) = \
                map_function(self, calculation)
            self.engine.children = [AMSJob(name=str(calculation.uid), molecule=plams_molecule, settings=plams_settings)]

        # If Simulation, create a loop of PLAMS jobs according to dependencies:
        elif list(self.dependencies.keys())[0] == "Simulation":
            job_list = []
            for calculation in self.dependencies["Simulation"]["Calculations"]:
                (plams_molecule, plams_settings) = map_function(self, calculation)
                job_list.append(AMSJob(name=str(calculation.uid), molecule=plams_molecule, settings=plams_settings))
            self.engine.children = job_list

    # OVERRIDE
    def _apply_updated(self, root_obj, buffer):
        """Updates the updated cuds in the engine."""
        # TODO: What should happen in the engine
        # when the user updates a certain cuds?
        # The given buffer contains all the updated CUDS object in a dictionary
#        map_results(self.engine, 'total_energy')

    # OVERRIDE
    def _apply_deleted(self, root_obj, buffer):
        """Deletes the deleted cuds from the engine."""
        # TODO: What should happen in the engine
        # when the user removes a certain cuds?
        # The given buffer contains all the deleted CUDS object in a dictionary

    @property
    def tarball(cls) -> str:
        return cls._tarball