"""Multiscale model for CO+Pt111"""
import tempfile
import warnings
from enum import Enum
from pathlib import Path
from typing import Union, List, TYPE_CHECKING, Optional
from uuid import UUID

from osp.core.namespaces import emmo, crystallography
from osp.core.session import CoreSession
from osp.core.utils import export_cuds
from osp.models.utils.general import (_download_file, _get_example_json,
                                      get_download, get_upload)
from osp.tools.io_functions import read_lattice, read_molecule
from pydantic import AnyUrl, Field, confloat, conint, conlist, root_validator
from pydantic.dataclasses import dataclass

from urllib.parse import quote, urlencode
from arcp import arcp_random

if TYPE_CHECKING:
    from osp.core.cuds import Cuds

STANDARD_XYZ = [
    ("4442d5c3-4b61-4b13-9bbb-fdf942776ca6", "CO_ads+Pt111.xyz"),
]


class ForceField(str, Enum):
    """Force field types in AMS."""

    OPLSAA = "OPLSAA"
    OPLS = "OPLS"
    AMBER = "AMBER"
    CHONSFPtClNi = "CHONSFPtClNi"
    ANAKINME = "ANAKINME"


class CalculationType(str, Enum):
    """Calculation types in AMS."""

    ProcessSearch = "ProcessSearch"
    WavefunctionOptimization = "WavefunctionOptimization"
    GeometryOptimization = "GeometryOptimization"
    PotentialEnergySurfaceScan = "PotentialEnergySurfaceScan"
    BindingSites = "BindingSites"
    StationaryPointCalculation = "StationaryPointCalculation"
    VibrationalFrequencyCalculation = "VibrationalFrequencyCalculation"


@dataclass
class PESExploration:
    """Data model for Potential Energy Surface (PES) Exploration."""

    force_field: ForceField = Field(
        ForceField.CHONSFPtClNi, description="Force field for the AMS calculation"
    )
    solver_type: str = Field("Direct", description="Type of solver used.")
    n_expeditions: conint(ge=1) = Field(
        30,
        description="""Sets the number of subsequent expeditions our job will consist of."""
        """Larger values result in a more comprehensive exploration of the potential energy surface, but will take more computational time.""",
    )
    n_explorers: conint(ge=1) = Field(
        3,
        description="""Sets the number of independent PES explorers dispatched as part of each expedition. 
                                      Larger values will result in a more comprehensive exploration of the potential energy surface, but will take more computational time.
                                       By default an appropriate number of explorers are executed in parallel.""",
    )
    max_energy: confloat(gt=0.0) = Field(2.0, description="Maximum energy in eV")
    max_distance: confloat(gt=0.0) = Field(
        3.8, description="Maximum distance cutoff from neighbors in Ångström"
    )
    random_seed: conint(ge=1) = Field(100, description="Random seed.")
    fixed_region: str = Field("surface", description="Fixed region of lattice.")
    reference_region: str = Field(
        "surface", description="Reference region in the lattice."
    )
    symmetry_check: str = Field(
        "T", description="Symmetry check for structure comparison."
    )
    molecule: Union[UUID, AnyUrl, Path] = Field(
        ...,
        description="""
        UUID of the cache-upload or url/system path to molecule input.""",
    )
    lattice: Union[UUID, AnyUrl, Path] = Field(
        ...,
        description="""
        UUID of the cache-upload or url/system path to lattica input.""",
    )

    def __post_init_post_parse__(self):
        with CoreSession() as session:
            self._make_model()
        file = tempfile.NamedTemporaryFile(suffix=".ttl", delete=False)
        export_cuds(session, file.name)
        self._file = file.name
        try:
            self._uuid = get_upload(file)
        except Exception as error:
            self._uuid = None
            message = (
                message
            ) = f"The graph of the model could not be stored at the minio-instance: {error.args}"
            warnings.warn(message)
        self._session = session

    def _make_model(self):
        calculation = emmo.ProcessSearch()

        if isinstance(self.molecule, UUID):
            xyz_file = get_download(str(self.molecule), as_file=True)
        elif isinstance(self.molecule, AnyUrl):
            xyz_file = _download_file(self.molecule, as_file=True)
        else:
            xyz_file = self.molecule
        molecule = read_molecule(xyz_file)

        if isinstance(self.lattice, UUID):
            xyz_file = get_download(str(self.lattice), as_file=True)
        elif isinstance(self.lattice, AnyUrl):
            xyz_file = _download_file(self.lattic, as_file=True)
        else:
            xyz_file = self.lattice
        lattice = read_lattice(xyz_file)

        forcefield = emmo[self.force_field.value]()

        solver = emmo.Solver()
        type_of_solver = emmo.Symbol(hasSymbolData=self.solver_type)
        solver.add(type_of_solver, rel=emmo.hasPart)

        fixed_region = emmo.FixedRegion(hasSymbolData=self.fixed_region)

        num_expeditions = emmo.NumberOfExpeditions(hasNumericalData=self.n_expeditions)
        num_explorers = emmo.NumberOfExplorers(hasNumericalData=self.n_explorers)

        max_energy = emmo.MaximumEnergy()
        energy_value = emmo.Real(hasNumericalData=self.max_energy)
        energy_unit = emmo.ElectronVolt(hasSymbolData="eV")
        max_energy.add(energy_value, rel=emmo.hasQuantityValue)
        max_energy.add(energy_unit, rel=emmo.hasReferenceUnit)

        max_distance = emmo.NeighborCutoff()
        distance_value = emmo.Real(hasNumericalData=self.max_distance)
        distance_unit = emmo.Ångström(hasSymbolData="Å")
        max_distance.add(distance_value, rel=emmo.hasQuantityValue)
        max_distance.add(distance_unit, rel=emmo.hasReferenceUnit)

        ref_region = emmo.ReferenceRegion(hasSymbolData=self.reference_region)
        random_seed = emmo.RandomSeed(hasNumericalData=self.random_seed)

        symmetry_check = emmo.CheckSymmetry(hasNumericalData="T")
        calculation.add(
            molecule,
            lattice,
            solver,
            fixed_region,
            num_expeditions,
            forcefield,
            num_explorers,
            max_energy,
            max_distance,
            ref_region,
            random_seed,
            symmetry_check,
            rel=emmo.hasInput,
        )
        return calculation

    @property
    def session(self) -> "CoreSession":
        return self._session

    @property
    def cuds(cls):
        return cls._session.load(cls._session.root).first()

    @property
    def uuid(cls):
        return cls._uuid

    @property
    def file(cls):
        return cls._file
    
    class Config:
        use_enum_values = True
        schema_extra = {
            "example": _get_example_json("pesexploration.json", STANDARD_XYZ)
        }

@dataclass
class BindingSite:
    """Data model for a binding site calculation based on a previous PESExploration"""

    n_expeditions: conint(ge=1) = Field(
        30,
        description="""Sets the number of subsequent expeditions our job will consist of."""
        """Larger values result in a more comprehensive exploration of the potential energy surface, but will take more computational time.""",
    )
    n_explorers: conint(ge=1) = Field(
        3,
        description="""Sets the number of independent PES explorers dispatched as part of each expedition. 
                                      Larger values will result in a more comprehensive exploration of the potential energy surface, but will take more computational time.
                                       By default an appropriate number of explorers are executed in parallel.""",
    )
    symmetry_check: str = Field(
        "F", description="Symmetry check for structure comparison."
    )


    def __post_init_post_parse__(self):
        with CoreSession() as session:
            self._make_model()
        self._session = session

    def _make_model(self):
        calculation = emmo.BindingSites()
        num_expeditions = emmo.NumberOfExpeditions(hasNumericalData=self.n_expeditions)
        num_explorers = emmo.NumberOfExplorers(hasNumericalData=self.n_explorers)
        symmetry_check = emmo.CheckSymmetry(hasNumericalData=self.symmetry_check)
        calculation.add(
            num_expeditions,
            num_explorers,
            symmetry_check,
            rel=emmo.hasInput,
        )
        return calculation

    @property
    def session(self) -> "CoreSession":
        return self._session

    @property
    def cuds(cls):
        return cls._session.load(cls._session.root).first()

@dataclass
class ZGBModel:
    """Data model for running ZGB model after PESExploration and BindingSite Calculation"""

    random_seed: int = Field(
        ..., description="""
        The integer seed of the random number generator.
        """
    )

    temperature: confloat(gt=0, allow_inf_nan=False) = Field(
        ..., description="""
        Temperature (K) under which the system is simulated.
        """
    )


    pressure: confloat(gt=0, allow_inf_nan=False) = Field(
        ..., description="""
        The pressure (bar) under which the system is simulated.
        """
    )

    max_time: confloat(ge=0, allow_inf_nan=True) = Field(
       ..., description="""
       The maximum allowed simulated time interval in seconds (time ranges from 0.0 to the maximum time in a
       simulation).
       """
    )

    n_gas_species: int = Field(
        ..., description="""
        The number of gas species in the chemistry.
        """
    )

    gas_specs_names: conlist(item_type=str) = Field(
        ..., description="""
        A list[str1, str2, ...] with the names of the gas species.
        There should be as many strings following the keyword as the number of gas species
        specified with keyword n_gas_species."
        """
    )

    gas_molar_fracs: conlist(item_type=confloat(ge=0, allow_inf_nan=False)) = Field(
        ..., description="""
        A list[str1, str2, ...] with the molar fractions (dim/less) of the gas species
        in the gas phase. There should be as many reals following this keyword as the number of gas
        species specified with keyword n_gas_species. The ordering of these values should be
        consistent with the order used in keyword gas_specs_names.
        """
    )

    snapshots: list = Field(
       ..., description="""
        Determines how often snapshots of the lattice state will be written to output file
        history_output.txt.
        """
    )

    species_numbers: list = Field(
       ..., description="""
        Determines how often information about the number of gas and surface species,
        as well as the energy of the current lattice configuration) will be written
        to specnum_output.txt
        """
    )

    process_statistics: list = Field(
       ..., description="""
        Determines how often statistical information about the occurrence of elementary
        events will be written to output file procstat_output.txt
        """
    )

    def __post_init_post_parse__(self):
        with CoreSession() as session:
            calculation = emmo.MesoscopicCalculation()
            calculation.add(*self._make_model(), rel=emmo.hasInput)
        self._session = session

    @root_validator
    def validate_all(cls, values):
        recording_options = ["off",
                             "on time",
                             "on event",
                             "on elemevent",
                             "on logtime",
                             "on realtime"]    
        if values["snapshots"][0] not in recording_options:
            raise ValueError(
                f"""Wrong recording option in 'snapshots', use one of the following:
                    {[recording_options]}.""")
        return values
    
    def _make_model(self) -> "List[Cuds]":
        return [
            *self._make_settings(),
            *self._make_gas_species(),
            self._make_snapshots(),
            self._make_species_numbers(),
            self._make_proc_stats()
        ]

    def _make_settings(self) -> "List[Cuds]":
        random_seed = emmo.RandomSeed(hasNumericalData=self.random_seed)

        temperature_float = emmo.Real(hasNumericalData=self.temperature)
        temperature_unit = emmo.Kelvin(hasSymbolData='K')
        temperature = emmo.ThermodynamicTemperature()
        temperature.add(temperature_unit, rel=emmo.hasReferenceUnit)
        temperature.add(temperature_float, rel=emmo.hasQuantityValue)

        pressure_float = emmo.Real(hasNumericalData=self.pressure)
        pressure_unit = emmo.Pascal(hasSymbolData='Pa')
        pressure = emmo.Pressure()
        pressure.add(pressure_unit, rel=emmo.hasReferenceUnit)
        pressure.add(pressure_float, rel=emmo.hasQuantityValue)

        time_float = emmo.Real(hasNumericalData=self.max_time)
        time_unit = emmo.Second(hasSymbolData='s')
        max_time = emmo.MaximumTime()
        max_time.add(time_unit, rel=emmo.hasReferenceUnit)
        max_time.add(time_float, rel=emmo.hasQuantityValue)

        return [random_seed, temperature, pressure, max_time]

    def _make_snapshots(self) -> "Cuds":
        recording_option = self.snapshots[0]
        iri = self._make_arcp("snapshots", query=dict(jsonpath=[["snapshot", str(0)]]))
        snap = emmo.Snapshots(hasSymbolData=recording_option, iri=iri)
        if recording_option != "off":

            if recording_option == "on logtime":
                cuds_object = emmo.Array()
                iri = self._make_arcp("snapshots", query=dict(jsonpath=[["snapshots", str(1)]]))
                time_float_1 = emmo.Real(hasNumericalData=self.snapshots[1], iri=iri)
                iri = self._make_arcp("snapshots", query=dict(jsonpath=[["snapshots", str(2)]]))
                time_float_2 = emmo.Real(hasNumericalData=self.snapshots[2], iri=iri)
                cuds_object.add(time_float_1, time_float_2, rel=emmo.hasSpatialPart)

            else:
                iri = self._make_arcp("snapshots", query=dict(jsonpath=[["snapshots", str(3)]]))
                cuds_object = emmo.Real(hasNumericalData=self.snapshots[1], iri=iri)

            snap.add(cuds_object, rel=emmo.hasSpatialPart)

        return snap

    def _make_species_numbers(self) -> "Cuds":
        recording_option = self.species_numbers[0]
        iri = self._make_arcp("species_numbers", query=dict(jsonpath=[["species_numbers", str(0)]]))
        numbers = emmo.SpeciesNumbers(hasSymbolData=recording_option)
        if recording_option != "off":

            if recording_option == "on logtime":
                cuds_object = emmo.Array()
                iri = self._make_arcp("species_numbers", query=dict(jsonpath=[["snapshots", str(1)]]))
                time_float_1 = emmo.Real(hasNumericalData=self.snapshots[1], iri=iri)
                iri = self._make_arcp("species_numbers", query=dict(jsonpath=[["snapshots", str(2)]]))
                time_float_2 = emmo.Real(hasNumericalData=self.snapshots[2], iri=iri)
                cuds_object.add(time_float_1, time_float_2, rel=emmo.hasSpatialPart)

            else:
                iri = self._make_arcp("species_numbers", query=dict(jsonpath=[["snapshots", str(3)]]))
                cuds_object = emmo.Real(hasNumericalData=self.snapshots[1], iri=iri)

            numbers.add(cuds_object, rel=emmo.hasSpatialPart)
        return numbers

    def _make_proc_stats(self) -> "Cuds":
        recording_option = self.process_statistics[0]
        iri = self._make_arcp("process_statistics",
                        query=dict(jsonpath=[["process_statistics", str(0)]]))
        stats = emmo.ProcessStatistics(hasSymbolData=recording_option, iri=iri)
        if recording_option != "off":

            if recording_option == "on logtime":
                cuds_object = emmo.Array()
                iri = self._make_arcp("process_statistics", query=dict(jsonpath=[["snapshots", str(1)]]))
                time_float_1 = emmo.Real(hasNumericalData=self.snapshots[1], iri=iri)
                iri = self._make_arcp("process_statistics", query=dict(jsonpath=[["snapshots", str(2)]]))
                time_float_2 = emmo.Real(hasNumericalData=self.snapshots[2], iri=iri)
                cuds_object.add(time_float_1, time_float_2, rel=emmo.hasSpatialPart)

            else:
                iri = self._make_arcp("process_statistics", query=dict(jsonpath=[["snapshots", str(3)]]))
                cuds_object = emmo.Real(hasNumericalData=self.snapshots[1], iri=iri)

            stats.add(cuds_object, rel=emmo.hasSpatialPart)

        return stats

    def _make_gas_species(self) -> "List[Cuds]":

        cuds = []
        for species in range(self.n_gas_species):

            gas_species = emmo.GasSpecies()
            molar_fraction = emmo.AmountFraction()
            symbol = emmo.ChemicalElement(hasSymbolData=self.gas_specs_names[species])

            iri = self._make_arcp("gas_molar_fracs",
                            query=dict(jsonpath=[["gas_molar_fracs", str(species)]]))
            molar_fraction_float = emmo.Real(
                hasNumericalData=self.gas_molar_fracs[species], iri=iri)
            molar_fraction.add(molar_fraction_float, rel=emmo.hasQuantityValue)

            gas_species.add(molar_fraction, symbol, rel=emmo.hasProperty)
            cuds.append(gas_species)

        return cuds

    def _make_arcp(self, *args, **kwargs):
        if kwargs.get("query"):
            query = kwargs.pop("query")
            for key, value in query.items():
                query[key] = [quote(".".join(item)) for item in value]
            query = urlencode(query, doseq=True)
        return arcp_random(*args, **kwargs, query=query)

    @property
    def session(self) -> "CoreSession":
        return self._session

    @property
    def cuds(cls):
        return cls._session.load(cls._session.root).first()


@dataclass
class COMolarFractionRange:
    """Range of molecular fractions of CO for adaptive design procedure."""

    min: confloat(ge=0.0, le=1.0) = Field(0.2, description="Minimum fraction of the named molecule.")
    max: confloat(ge=0.0, le=1.0) = Field(0.8, description="Maximum fraction of the named molecule")
    num: int = Field(5, description="Number of elements in the range of molar fractions.")

    def __post_init_post_parse__(self):
        with CoreSession() as session:
            molecule = emmo.MolecularGeometry()
            name = emmo.ChemicalName(hasSymbolData="CO")
            vector = emmo.Vector()
            frac = emmo.AmountConcentration()
            maximum = emmo.Real(hasNumericalData=self.max)
            minimum = emmo.Real(hasNumericalData=self.min)
            length = emmo.Real(hasNumericalData=self.num)
            vector.add(maximum, rel=emmo.hasMaximumValue)
            vector.add(minimum, rel=emmo.hasMinimumValue)
            vector.add(length, rel=emmo.hasVectorLength)
            frac.add(vector, rel=emmo.hasSign)
            molecule.add(frac, rel=emmo.hasQuantitativeProperty)
            molecule.add(name, rel=emmo.hasProperty)
        file = tempfile.NamedTemporaryFile(suffix=".ttl", delete=False)
        export_cuds(session, file.name)
        self._file = file.name
        try:
            self._uuid = get_upload(file)
        except Exception as error:
            self._uuid = None
            message = (
                message
            ) = f"The graph of the model could not be stored at the minio-instance: {error.args}"
            warnings.warn(message)
        self._session = session

    @property
    def session(self) -> "CoreSession":
        return self._session

    @property
    def cuds(cls):
        return cls._session.load(cls._session.root).first()

    @property
    def uuid(cls):
        return cls._uuid

    @property
    def file(cls):
        return cls._file



@dataclass
class COPt111MesoscaleModel:
    """General Model for multiscale simulation for PESExploration, 
    Binding Site Calculation and ZGB Model"""

    pes_exploration: PESExploration = Field(
        ..., description="AMS data model for PESExploration."
    )
    binding_site: BindingSite = Field(
        ...,
        description="""data model for binding site calculation
        based on the previous PESExploraion."""
    )
    zgb_model: ZGBModel = Field(
        ..., description="ZGB model for mesoscopic scale."
    )

    adp: Optional[COMolarFractionRange] = Field(
        None,
        description="""Molar fractions of CO
        for the adaptive design procedure"""
        )

    def __post_init_post_parse__(self):
        with CoreSession() as session:
            workflow = emmo.Workflow()
            self.pes_exploration.cuds.add(self.binding_site.cuds, rel=emmo.hasSpatialNext)
            self.binding_site.cuds.add(self.zgb_model.cuds, rel=emmo.hasSpatialNext)
            for oclass in [
                    emmo.ForceFieldIdentifierString,
                    emmo.Solver,
                    emmo.FixedRegion,
                    emmo.MaximumEnergy,
                    emmo.NeighborCutoff,
                    emmo.ReferenceRegion,
                    emmo.RandomSeed,
                    emmo.MolecularGeometry,
                    crystallography.UnitCell
                ]:
                input_cuds = self.pes_exploration.cuds.get(oclass=oclass, rel=emmo.hasInput)
                self.binding_site.cuds.add(input_cuds.pop(), rel=emmo.hasInput)
            workflow.add(self.pes_exploration.cuds, rel=emmo.hasSpatialFirst)
            workflow.add(self.binding_site.cuds, rel=emmo.hasSpatialDirectPart)
            if self.adp:
                apd = emmo.AdaptiveDesignProcedure()
                apd.add(self.adp.cuds, rel=emmo.hasInput)
                self.zgb_model.cuds.add(apd, rel=emmo.hasSpatialNext)
                workflow.add(apd, rel=emmo.hasSpatialLast)
                workflow.add(self.zgb_model.cuds, rel=emmo.hasSpatialDirectPart)
            else:
                workflow.add(self.zgb_model.cuds, rel=emmo.hasSpatialLast)
        file = tempfile.NamedTemporaryFile(suffix=".ttl", delete=False)
        export_cuds(session, file.name)
        self._file = file.name
        try:
            self._uuid = get_upload(file)
        except Exception as error:
            self._uuid = None
            message = (
                message
            ) = f"The graph of the model could not be stored at the minio-instance: {error.args}"
            warnings.warn(message)
        self._session = session

    @property
    def session(self) -> "CoreSession":
        return self._session

    @property
    def cuds(cls):
        return cls._session.load(cls._session.root).first()

    @property
    def uuid(cls):
        return cls._uuid

    @property
    def file(cls):
        return cls._file

    class Config:
        """Pydantic Config"""
    
        schema_extra = {
            "example": _get_example_json("co_pt111_meso.json", STANDARD_XYZ)
        }

