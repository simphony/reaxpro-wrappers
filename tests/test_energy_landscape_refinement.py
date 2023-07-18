"""Pytest for the pydantic model of the energy landscape refinement"""
import os
from typing import TYPE_CHECKING

PATH = os.path.dirname(__file__)

if TYPE_CHECKING:
    from osp.models.ams.energy_landscape_refinement import EnergyLandscapeRefinement


def _check_model(
    model: "EnergyLandscapeRefinement",
    n_molecules: int = 5,
    n_pathways: int = 5,
    n_transitionstates: int = 2,
) -> None:
    """Helper function for checking the consistency of the CUDS"""
    from osp.core.namespaces import emmo

    assert model.cuds.is_a(emmo.LandscapeRefinement)
    landscape = model.cuds.get(oclass=emmo.EnergyLandscape)
    assert len(landscape) == 1
    landscape = landscape.pop()
    pathways = landscape.get(oclass=emmo.ReactionPathway)
    assert len(pathways) == n_pathways
    react_pathways = [
        pathway
        for pathway in pathways
        if pathway.get(oclass=emmo.ChemicalReactionEquationProduct)
        and pathway.get(oclass=emmo.ChemicalReactionEquationReactant)
    ]
    mole_pathways = [
        pathway
        for pathway in pathways
        if pathway.get(oclass=emmo.MolecularGeometry)
        and not pathway.get(oclass=emmo.TransitionStateGeometry)
    ]
    trans_pathways = [
        pathway
        for pathway in pathways
        if pathway.get(oclass=emmo.TransitionStateGeometry)
    ]
    assert len(react_pathways) == n_transitionstates
    assert len(mole_pathways) == n_molecules - n_transitionstates
    assert len(trans_pathways) == n_transitionstates
    molecules = set()
    for pathway in react_pathways:
        reactant = pathway.get(oclass=emmo.ChemicalReactionEquationReactant).pop()
        product = pathway.get(oclass=emmo.ChemicalReactionEquationProduct).pop()
        reactant = reactant.get(oclass=emmo.MolecularGeometry)
        product = product.get(oclass=emmo.MolecularGeometry)
        assert len(reactant) == 1
        assert len(product) == 1
        molecules.update({reactant.pop(), product.pop()})
    for pathway in mole_pathways:
        molecule = pathway.get(oclass=emmo.MolecularGeometry)
        assert len(molecule) == 1
        molecules.update({molecule.pop()})
    for pathway in trans_pathways:
        transition = pathway.get(oclass=emmo.TransitionStateGeometry)
        assert len(transition) == 1
    assert len(molecules) == n_molecules - n_transitionstates


def test_energy_landscape_models():
    """Load the energy landscape from an filesyste-path."""

    from osp.models.ams.energy_landscape_refinement import EnergyLandscapeRefinement

    from .conftest import store_file

    uuids = []
    for i in range(1, 6):
        state = os.path.join(PATH, "test_files", f"state{i}.xyz")
        with open(state, "rb") as file:
            key = store_file(file)
        uuids.append(key)

    config = {
        "pathways": [
            {
                "reactant": {"xyz_file": uuids[0], "charge": -1},
                "transition": {"xyz_file": uuids[2], "charge": -1},
                "product": {"xyz_file": uuids[1], "charge": -1},
            },
            {
                "reactant": {"xyz_file": uuids[1], "charge": -1},
                "transition": {"xyz_file": uuids[4], "charge": -1},
                "product": {"xyz_file": uuids[3], "charge": -1},
            },
        ]
    }
    model = EnergyLandscapeRefinement(**config)
    _check_model(model)


def test_energy_landscape_model_filepath():
    from osp.models.ams.energy_landscape_refinement import EnergyLandscapeRefinement

    state1 = os.path.join(PATH, "test_files", f"state1.xyz")
    state2 = os.path.join(PATH, "test_files", f"state2.xyz")
    state3 = os.path.join(PATH, "test_files", f"state3.xyz")
    state4 = os.path.join(PATH, "test_files", f"state4.xyz")
    state5 = os.path.join(PATH, "test_files", f"state5.xyz")
    config = {
        "pathways": [
            {
                "reactant": {"xyz_file": state1, "charge": -1},
                "transition": {"xyz_file": state3, "charge": -1},
                "product": {"xyz_file": state2, "charge": -1},
            },
            {
                "reactant": {"xyz_file": state2, "charge": -1},
                "transition": {"xyz_file": state5, "charge": -1},
                "product": {"xyz_file": state4, "charge": -1},
            },
        ]
    }
    model = EnergyLandscapeRefinement(**config)
    _check_model(model)


def test_energy_landscape_model_filepath_wo_transition():
    from osp.models.ams.energy_landscape_refinement import EnergyLandscapeRefinement

    state1 = os.path.join(PATH, "test_files", f"state1.xyz")
    state2 = os.path.join(PATH, "test_files", f"state2.xyz")
    state4 = os.path.join(PATH, "test_files", f"state4.xyz")
    state5 = os.path.join(PATH, "test_files", f"state5.xyz")
    config = {
        "pathways": [
            {
                "reactant": {"xyz_file": state1, "charge": -1},
                "product": {"xyz_file": state2, "charge": -1},
            },
            {
                "reactant": {"xyz_file": state2, "charge": -1},
                "transition": {"xyz_file": state5, "charge": -1},
                "product": {"xyz_file": state4, "charge": -1},
            },
        ]
    }
    model = EnergyLandscapeRefinement(**config)
    _check_model(model, n_molecules=4, n_pathways=4, n_transitionstates=1)


def test_energy_landscape_model_example():
    from osp.models.ams.energy_landscape_refinement import EnergyLandscapeRefinement
    from osp.models.utils.general import _get_example

    example = _get_example("ams", "energy_landscape.json")
    assert isinstance(example, dict)
    assert example != {}
    assert example == EnergyLandscapeRefinement.Config.schema_extra["example"]

    model = EnergyLandscapeRefinement(**example)
    _check_model(model)
