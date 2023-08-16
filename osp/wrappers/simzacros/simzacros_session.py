"""Describe Zacros Session class."""
from osp.core.namespaces import emmo, crystallography
from osp.core.session import SimWrapperSession
from osp.tools.graph_functions import graph_wrapper_dependencies
from osp.core.cuds import Cuds
from osp.tools.io_functions import raise_error
from osp.tools.mapping_functions import map_function, map_results
from osp.core.utils import simple_search as search
from osp.models.utils.general import get_download, get_upload
from scm.plams import config
from uuid import uuid4, UUID
import logging
import scm.pyzacros as pz
import numpy as np
import os
# from osp.core.utils import pretty_print

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class SimzacrosSession(SimWrapperSession):
    """Describe Zacros Session class."""

    def __init__(self, engine=None, **kwargs):
        """Initialise SimamsSession."""
        if engine is None:
            pz.init()
            config.job.runscript.nproc = int(os.environ.get("REAXPRO_N_PROCESSES")) or 1
            self.workdir = config.get("default_jobmanager").workdir
            self.jobname = str(uuid4())
            self.engine = "Zacros"
            self.calculation = None
            self.mechanism = None
            self.cluster = None
            self.adp = []

        super().__init__(engine)

    def __str__(self):
        """To overwrite the private str method. Not advised, but here it is."""
        # TODO: Define the output of str(SomeSimulationSession())
        return "pyZacros Wrapper Session"

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
        if self.adp:
            import adaptiveDesignProcedure as adp

            def get_rate( conditions ):

                logger.info(f"Requesting conditions for x_CO: {conditions}")
                #---------------------------------------
                # Zacros calculation
                #---------------------------------------

                ps_params = pz.ZacrosParametersScanJob.Parameters()
                ps_params.add( 'x_CO', 'molar_fraction.CO', [ cond[0] for cond in conditions ] )
                ps_params.add( 'x_O2', 'molar_fraction.O2', lambda params: 1.0-params['x_CO'] )

                ps_job = pz.ZacrosParametersScanJob( reference=pz_job, parameters=ps_params )

                #---------------------------------------
                # Running the calculations
                #---------------------------------------
                results = ps_job.run()

                if not results.job.ok():
                    logger.error(
                        f"""Something went wrong requesting
                        parameters from {pz_job}"""
                    )

                #---------------------------------------
                # Collecting the results
                #---------------------------------------
                data = np.nan*np.empty((len(conditions),1))
                results_dict = results.turnover_frequency()

                for i in range(len(results_dict)):
                    data[i,0] = results_dict[i]['turnover_frequency']['CO2']

                return data
            output_var = [{'name':'TOF_CO2'}]
            adp_path = os.path.join(self.workdir,'adp.results')
            adpML = adp.adaptiveDesignProcedure(
                self.adp, 
                output_var, 
                get_rate,
                outputDir=adp_path,
                randomState=10
            )
            adpML.createTrainingDataAndML()

            uuid = get_upload(adpML.forestFileForCFD)
            pkl = emmo.PKLFile(uid=UUID(uuid))
            self.adp_cuds.add(pkl, rel=emmo.hasOutput)


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

        dependencies = graph_wrapper_dependencies(root_obj, oclass=emmo.MesoscopicCalculation)
        if "Simulation" in dependencies:
            calculations = dependencies["Simulation"]["Calculations"]
        elif "Calculation" in dependencies:
            calculations = dependencies["Calculation"]
        else:
            calculations = []

        for calculation in calculations:
            if calculation.is_a(emmo.MesoscopicCalculation):
                if self.calculation:
                    raise_error(file=os.path.basename(__file__),
                        function="_apply_added method",
                        type='NameError',
                        message='More than one emmo.MesoscopicCalculation defined'
                        ' in the Wrapper object.')
                self._process_calculation(calculation)
            elif calculation.is_a(emmo.AdaptiveDesignProcedure):
                if self.adp:
                    raise_error(file=os.path.basename(__file__),
                        function="_apply_added method",
                        type='NameError',
                        message='More than one emmo.AdaptiveDesignProcedure defined'
                        ' in the Wrapper object.')        
                self._process_adp(calculation)

    def _process_calculation(self, calculation: Cuds) -> None:
        # Calculation
        self.calculation = calculation

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
        if len(search_lattice) > 1:
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

        # Cluster
        search_cluster = \
            search.find_cuds_objects_by_oclass(emmo.ClusterExpansion,
                                               self.calculation, emmo.hasInput)
        if search_cluster:
            self.cluster = search_cluster

    def _process_adp(self, calculation: Cuds) -> None:
        self.adp_cuds = calculation
        molecule = calculation.get(oclass=emmo.MolecularGeometry, rel=emmo.hasInput)
        if not molecule:
            raise ValueError(f"No <{emmo.MolecularGeometry}> found in <{calculation}>.")
        elif len(molecule) > 1:
            raise ValueError(f"More than 1 <{emmo.MolecularGeometry}> found in <{calculation}>.")
        molecule = molecule.pop()
        name = molecule.get(oclass=emmo.ChemicalName, rel=emmo.hasProperty)
        if not name:
            raise ValueError(f"No <{emmo.ChemicalName}> found in <{molecule}>.")
        elif len(name) > 1:
            raise ValueError(f"More than 1 <{emmo.ChemicalName}> found in <{molecule}>.")
        name = name.pop()
        if not name.hasSymbolData == "CO":
            raise ValueError(f"<{molecule}> does not have `CO` as <{name}> via <{emmo.hasSymbolData}>.")
        frac = molecule.get(oclass=emmo.AmountConcentration, rel=emmo.hasQuantitativeProperty)
        if not frac:
            raise ValueError(f"<{molecule}> does not have any <{emmo.AmountConcentration}>.")
        elif len(frac) > 1:
            raise ValueError(f"<{molecule}> has more than 1 <{emmo.AmountConcentration}>.")
        frac = frac.pop()
        vec = frac.get(oclass=emmo.Vector, rel=emmo.hasSign)
        if not vec:
            raise ValueError(f"<{frac}> does not have any <{emmo.Vector}>.")
        elif len(vec) > 1:
            raise ValueError(f"<{frac}> has more than 1 <{emmo.Vector}>.")
        vec = vec.pop()
        maximum = vec.get(oclass=emmo.Real, rel=emmo.hasMaximumValue)
        minimum = vec.get(oclass=emmo.Real, rel=emmo.hasMinimumValue)
        length = vec.get(oclass=emmo.Real, rel=emmo.hasVectorLength)
        if not minimum or not maximum or not length:
            raise ValueError(
                f"""Any of the following quantities are empty: 
                Minimum {minimum}, maximum: {maximum}, vector length {length}"""
                )
        elif len(minimum) > 1 or len(maximum) > 1 or len(length) > 1:
            raise ValueError(
                f"""Any of the following quantities has more then 1 times:
                Minimum {minimum}, maximum: {maximum}, vector length {length}"""
                )
        maximum, minimum, length = maximum.pop(), minimum.pop(), length.pop()
        self.adp.append(
            {
                "name": f"x_{name.hasSymbolData}",
                "min": float(minimum.hasNumericalData),
                "max": float(maximum.hasNumericalData),
                "num": int(length.hasNumericalData)
            }
        )


    # OVERRIDE
    def _apply_updated(self, root_obj, buffer):
        """Updates the updated cuds in the engine."""
        # TODO: What should happen in the engine
        # when the user updates a certain cuds?
        # The given buffer contains all the updated CUDS object in a dictionary
        # map_results(self.engine, 'total_energy')

    # OVERRIDE
    def _apply_deleted(self, root_obj, buffer):
        """Deletes the deleted cuds from the engine."""
        # TODO: What should happen in the engine
        # when the user removes a certain cuds?
        # The given buffer contains all the deleted CUDS object in a dictionary

    @property
    def tarball(cls) -> str:
        return cls._tarball