import tempfile
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional, Union
from uuid import UUID

from pydantic import AnyUrl, BaseModel, Field, root_validator
from pydantic.dataclasses import dataclass

from osp.core.namespaces import emmo
from osp.core.session import CoreSession
from osp.core.utils import export_cuds
from osp.models.utils.general import (_download_file, _get_example_json,
                                      get_download, get_upload)
from osp.tools.io_functions import GeometryType, read_molecule

STANDARD_XYZ = [
    ("b66ff7f5-979f-4e38-8ad3-e6bce42145d5", "state1.xyz"),
    ("76dd6a5c-0e80-4b82-af2e-0b2f884af1b0", "state3.xyz"),
    ("d1c1b58b-50ef-4e1b-9776-6ad90b20fe98", "state2.xyz"),
    ("8944c8d3-e7a4-453d-a179-d9126d1c5d61", "state5.xyz"),
    ("934651a0-badd-43a3-98fd-ed67265b9669", "state4.xyz"),
]


if TYPE_CHECKING:
    from typing import Dict

    from osp.core.cuds import Cuds


class Molecule(BaseModel):
    xyz_file: Union[UUID, AnyUrl, Path] = Field(
        ...,
        description="UUID of the cache-upload or url/system path to xyz file for the pathway of the molecule.",
    )
    charge: int = Field(..., description="Electric charge of the moluecole in eV.")


class MoleculeReaction(BaseModel):
    reactant: Molecule = Field(
        ..., description="Reactant molecule within the chemical reaction."
    )
    transition: Optional[Molecule] = Field(
        None,
        description="Transition state geometry  of the molecule within the chemical reaction.",
    )
    product: Molecule = Field(
        ..., description="Product molecule within the chemical reaction."
    )

    @root_validator(allow_reuse=True)
    def validate_all(cls, values):
        if not values.get("transition"):
            if values["reactant"].xyz_file == values["product"].xyz_file:
                raise ValueError("Reactant and product are the same file references!")
        else:
            for file1, file2 in (
                ("reactant", "transition"),
                ("reactant", "product"),
                ("transition", "product"),
            ):
                if values[file1].xyz_file == values[file2].xyz_file:
                    raise ValueError(
                        file1 + " and " + file2 + "are the same file references!"
                    )
        return values


@dataclass
class EnergyLandscapeRefinement:
    """Run AMS LandscapeRefinment calculation."""

    pathways: List[MoleculeReaction] = Field(
        ..., description="""List of molecule pathways included for the refinement."""
    )

    def __post_init_post_parse__(self):
        with CoreSession() as session:
            calculation = emmo.LandscapeRefinement()
            model = emmo.DFTB()
            calculation.add(self._make_landscape(), model, rel=emmo.hasInput)
        file = tempfile.NamedTemporaryFile(suffix=".ttl",delete=False)
        export_cuds(session, file.name)
        self._file = file.name
        try:
            self._uuid = get_upload(file)
        except Exception as error:
            self._uuid = None
            message = f"The graph of the model could not be stored at the minio-instance: {error.args}"
            warnings.warn(message)
        self._session = session

    def _make_landscape(self) -> "Cuds":
        energy_landscape = emmo.EnergyLandscape()
        pathways = self._make_reaction_paths()
        for index in range(len(pathways)):
            if index == 0:
                energy_landscape.add(pathways[index], rel=emmo.hasSpatialFirst)
            else:
                pathways[index-1].add(pathways[index], rel=emmo.hasSpatialNext)
                energy_landscape.add(pathways[index], rel=emmo.hasSpatialDirectPart)
            if index == len(pathways) + 1:
                energy_landscape.add(pathways[index], rel=emmo.hasSpatialLast)
        return energy_landscape

    def _make_reaction_paths(self) -> "List[Cuds]":
        reactions = self._make_reaction()
        paths = []
        for index in range(len(reactions)):
            if index == 0:
                molecule1 = reactions[index]["reactant"]
                path1 = emmo.ReactionPathway()
                path1.add(molecule1, rel=emmo.hasPart)
                paths.append(path1)
            else:
                molecule1 = reactions[index - 1]["product"]
            molecule2 = reactions[index].get("transition")
            molecule3 = reactions[index]["product"]
            reactant = emmo.ChemicalReactionEquationReactant()
            reactant.add(molecule1, rel=emmo.hasSpatialDirectPart)
            product = emmo.ChemicalReactionEquationProduct()
            product.add(molecule3, rel=emmo.hasSpatialDirectPart)

            path2 = emmo.ReactionPathway()
            path2.add(molecule3, rel=emmo.hasPart)
            paths.append(path2)

            if molecule2:
                path3 = emmo.ReactionPathway()
                path3.add(reactant, product, molecule2, rel=emmo.hasSpatialDirectPart)
                paths.append(path3)

        return paths

    def _make_reaction(self) -> "List[Dict[str, Cuds]]":
        init_pathway = self.pathways[0]
        cuds_pathways = [
            {
                "reactant": self._make_molecule(init_pathway.reactant),
                "product": self._make_molecule(init_pathway.product),
            }
        ]
        if init_pathway.transition:
            cuds_pathways[0]["transition"] = self._make_molecule(
                init_pathway.transition, TS=True
            )
        for index in range(1, len(self.pathways)):
            ipathway = self.pathways[index]
            reactant = cuds_pathways[-1]["product"]
            cuds_pathways.append(
                {
                    "product": self._make_molecule(ipathway.reactant),
                    "reactant": reactant,
                }
            )
            if ipathway.transition:
                cuds_pathways[-1]["transition"] = self._make_molecule(
                    ipathway.transition, TS=True
                )
        return cuds_pathways

    def _make_molecule(self, molecule: Molecule, TS: bool = False) -> "Cuds":
        if isinstance(molecule.xyz_file, UUID):
            xyz_file = get_download(str(molecule.xyz_file), as_file=True)
        elif isinstance(molecule.xyz_file, AnyUrl):
            xyz_file = _download_file(molecule.xyz_file, as_file=True)
        else:
            xyz_file = molecule.xyz_file
        cuds = read_molecule(xyz_file, geometry_type=GeometryType.XYZ, TS=TS)
        charge = emmo.ElectricCharge()
        integer = emmo.Integer(hasNumericalData=molecule.charge)
        charge.add(integer, rel=emmo.hasQuantityValue)
        cuds.add(charge, rel=emmo.hasProperty)
        return cuds

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
        """Pydantic Config for EnergyLandscapeRefinement"""

        schema_extra = {
            "example": _get_example_json("energy_landscape.json", STANDARD_XYZ)
        }
