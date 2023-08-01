"""Multiscale model for CO+Pt111"""
import tempfile
import warnings
from enum import Enum
from pathlib import Path
from typing import Union
from uuid import UUID

from osp.core.namespaces import emmo
from osp.core.session import CoreSession
from osp.core.utils import export_cuds
from osp.models.utils.general import (_download_file, _get_example_json,
                                      get_download, get_upload)
from osp.models.zacros.co_pyzacros import COpyZacrosModel
from osp.tools.io_functions import read_lattice, read_molecule
from pydantic import AnyUrl, BaseModel, Field, confloat, conint
from pydantic.dataclasses import dataclass


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
        ForceField.CHONSFPtClNi.value, description="Force field for the AMS calculation"
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

        forcefield = emmo[self.force_field]()

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
        num_expeditions = emmo.NumberOfExpeditions(hasNumericalData=self.n_expedition)
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
class COPt111MesoscaleModel(BaseModel):
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
    zgb_model: COpyZacrosModel = Field(
        ..., description="ZGB model for mesoscopic scale."
    )

    def __post_init_post_parse__(self):
        with CoreSession() as session:
            workflow = emmo.Workflow()
            workflow.add(self.pes_exploration.cuds, rel=emmo.hasSpatialFirst)
            self.pes_exploration.cuds.add(self.binding_site.cuds, rel=emmo.hasSpatialNext)
            workflow.add(self.binding_site.cuds, rel=emmo.hasSpatialDirectPart)
            self.binding_site.cuds.add(self.zacros_model.cuds, rel=emmo.hasSpatialNext)
            workflow.add(self.zgb_model.cuds, rel=emmo.hasSpatialLast)
            for oclass in [
                    "ForceFieldIdentifierString",
                    "Solver",
                    "FixedRegion",
                    "MaximumEnergy",
                    "NeighborCutoff",
                    "ReferenceRegion",
                    "RandomSeed"
                ]:
                input_cuds = self.pes_exploration.cuds.get(oclass=emmo[oclass], rel=emmo.hasInput)
                self.binding_site.cuds.add(input_cuds.pop(), rel=emmo.hasInput)
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
