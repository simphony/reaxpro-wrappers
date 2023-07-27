"""Describe Zacros Session class."""
from osp.core.namespaces import emmo, crystallography, cuba
from osp.core.session import SimWrapperSession
from osp.core.cuds import Cuds
from osp.tools.io_functions import raise_error
from osp.tools.mapping_functions import map_function, map_results
from osp.core.utils import simple_search as search
from osp.models.utils.general import get_download
import scm.pyzacros as pz
import os
# from osp.core.utils import pretty_print


class SimzacrosSession(SimWrapperSession):
    """Describe Zacros Session class."""

    def __init__(self, engine=None, **kwargs):
        """Initialise SimamsSession."""
        if engine is None:
            pz.init()
            self.engine = 'Zacros'
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

        (pz_settings, pz_lattice, pz_mechanism, pz_cluster_expansion) = \
            map_function(self, root_cuds_object, self.engine)
        pz_job = pz.ZacrosJob(settings=pz_settings, lattice=pz_lattice,
                              mechanism=pz_mechanism, cluster_expansion=pz_cluster_expansion)
        results = pz_job.run()
        self._tarball = map_results(pz_job, root_cuds_object)

        # Attributes to easily access (syntactic) info from results.
        self.get_reaction_network = results.get_reaction_network()
        self.provided_quantities_names = results.provided_quantities_names()
        self.provided_quantities = results.provided_quantities()
        self.number_of_lattice_sites = results.number_of_lattice_sites()
        self.gas_species_names = results.gas_species_names()
        self.surface_species_names = results.surface_species_names()
        self.site_type_names = results.site_type_names()
        self.number_of_snapshots = results.number_of_snapshots()
        self.number_of_process_statistics = results.number_of_process_statistics()
        self.elementary_steps_names = results.elementary_steps_names()

    # OVERRIDE
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

        self.input_path = '../tests/test_files'

        # Calculation
        search_calculation = \
            search.find_cuds_objects_by_oclass(
                           emmo.MesoscopicCalculation,
                           root_obj, cuba.relationship)
        #  Raise error if two calculations, for some reason, are found:
        if len(search_calculation) > 1:
            raise_error(file=os.path.basename(__file__),
                        function="_apply_added method",
                        type='NameError',
                        message='More than one emmo.MesoscopicCalculation defined'
                        ' in the Wrapper object.')
        if search_calculation:
            self.calculation = search_calculation[0]

        # Mechanism
        search_mechanism = \
            search.find_cuds_objects_by_oclass(
                           emmo.ChemicalReactionMechanism,
                           self.calculation, emmo.hasInput)
        #  Raise error if two mechanism, for some reason, are found:
        if len(search_mechanism) > 1:
            raise_error(file=os.path.basename(__file__),
                        function="_apply_added method",
                        type='NameError',
                        message='More than one emmo.ChemicalReactionMechanism defined'
                        ' in the Wrapper object.')
        if search_mechanism:
            self.mechanism = search_mechanism[0]

        # Lattice
        search_lattice = \
            search.find_cuds_objects_by_oclass(
                           crystallography.UnitCell,
                           self.calculation, emmo.hasInput)

        #  Raise error if two calculations, for some reason, are found:
        if len(search_calculation) > 1:
            raise_error(file=os.path.basename(__file__),
                        function="_apply_added method",
                        type='NameError',
                        message='More than one crystallography.UnitCell defined'
                        ' in the Wrapper object.')
        if search_lattice:
            lattice = search_lattice.pop()
            if "file://" in str(lattice.iri):
                split = str(lattice.iri).split("file://")
                self.lattice = split[-1]
            else:
                self.lattice = get_download(str(lattice.uid), as_file=True)
        else:
            # lattice loaded from semantic-based script:
            self.lattice = search_lattice[0].iri[49:]

        # Cluster
        search_cluster = \
            search.find_cuds_objects_by_oclass(emmo.ClusterExpansion,
                                               self.calculation, emmo.hasInput)
        if search_cluster:
            self.cluster = search_cluster

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