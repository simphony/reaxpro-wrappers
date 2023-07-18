"""Run AMS LandscapeRefinment calculation using ReaxPro ontology.""" 

from osp.core.namespaces import emmo, cuba
from osp.tools.io_functions import read_molecule
from osp.wrappers.simams.simams_session import SimamsSession
from osp.core.utils import pretty_print
from osp.core.utils import Cuds2dot
from osp.core.utils import simple_search as search
from osp.core.utils import import_cuds, export_cuds

calculation = emmo.LandscapeRefinement()
model = emmo.DFTB()

molecule1 = read_molecule('XYZ/state1.xyz')
molecule1_charge = emmo.ElectricCharge()
molecule1_integer = emmo.Integer(hasNumericalData=-1)
molecule1_charge.add(molecule1_integer, rel=emmo.hasQuantityValue)
molecule1.add(molecule1_charge, rel=emmo.hasProperty)

molecule2 = read_molecule('XYZ/state2.xyz')
molecule2_charge = emmo.ElectricCharge()
molecule2_integer = emmo.Integer(hasNumericalData=-1)
molecule2_charge.add(molecule2_integer, rel=emmo.hasQuantityValue)
molecule2.add(molecule2_charge, rel=emmo.hasProperty)

molecule3 = read_molecule('XYZ/state3.xyz')
molecule3_charge = emmo.ElectricCharge()
molecule3_integer = emmo.Integer(hasNumericalData=-1)
molecule3_charge.add(molecule3_integer, rel=emmo.hasQuantityValue)
molecule3.add(molecule3_charge, rel=emmo.hasProperty)

molecule4 = read_molecule('XYZ/state4.xyz')
molecule4_charge = emmo.ElectricCharge()
molecule4_integer = emmo.Integer(hasNumericalData=-1)
molecule4_charge.add(molecule4_integer, rel=emmo.hasQuantityValue)
molecule4.add(molecule4_charge, rel=emmo.hasProperty)

molecule5 = read_molecule('XYZ/state5.xyz')
molecule5_charge = emmo.ElectricCharge()
molecule5_integer = emmo.Integer(hasNumericalData=-1)
molecule5_charge.add(molecule5_integer, rel=emmo.hasQuantityValue)
molecule5.add(molecule5_charge, rel=emmo.hasProperty)

energy_landscape = emmo.EnergyLandscape()

path_1 = emmo.ReactionPathway()
path_1.add(molecule1, rel=emmo.hasPart)

path_2 = emmo.ReactionPathway()
path_2.add(molecule2, rel=emmo.hasPart)

path_3 = emmo.ReactionPathway() 
path_3_reactant = emmo.ChemicalReactionEquationReactant()
path_3_reactant.add(molecule1, rel=emmo.hasSpatialDirectPart)

path_3_TS = read_molecule('XYZ/state3.xyz', TS=True)
path_3_TS_charge = emmo.ElectricCharge()
path_3_TS_integer = emmo.Integer(hasNumericalData=-1)
path_3_TS_charge.add(molecule3_integer, rel=emmo.hasQuantityValue)
path_3_TS.add(molecule3_charge, rel=emmo.hasProperty)

path_3_product = emmo.ChemicalReactionEquationProduct()
path_3_product.add(molecule2, rel=emmo.hasSpatialDirectPart)
path_3.add(path_3_reactant, path_3_product, path_3_TS, rel=emmo.hasPart)
path_3.add(path_3, rel=emmo.hasPart)

path_4 = emmo.ReactionPathway()
path_4.add(molecule4, rel=emmo.hasPart)

path_5 = emmo.ReactionPathway()
path_5_reactant = emmo.ChemicalReactionEquationReactant()
path_5_reactant.add(molecule2, rel=emmo.hasSpatialDirectPart)
path_5_product = emmo.ChemicalReactionEquationProduct()
path_5_product.add(molecule4, rel=emmo.hasSpatialDirectPart)

path_5_TS = read_molecule('XYZ/state5.xyz')
path_5_TS_charge = emmo.ElectricCharge()
path_5_TS_integer = emmo.Integer(hasNumericalData=-1)
path_5_TS_charge.add(molecule5_integer, rel=emmo.hasQuantityValue)
path_5_TS.add(molecule5_charge, rel=emmo.hasProperty)

path_5.add(path_5_reactant, path_5_product, path_5_TS, rel=emmo.hasPart)
path_5.add(path_5, rel=emmo.hasPart)

energy_landscape.add(path_1, path_2, path_3, path_4, path_5, rel=emmo.hasPart)
calculation.add(energy_landscape, model, rel=emmo.hasInput)
pretty_print(energy_landscape)


with SimamsSession() as sess:
    reaxpro_wrapper = cuba.Wrapper(session=sess)
    reaxpro_wrapper.add(calculation,
                        rel=emmo.hasPart)
    reaxpro_wrapper.session.run()
#
## Post-processing:
Cuds2dot(calculation).render()
search_calculation = \
        search.find_cuds_objects_by_oclass(
                           emmo.Calculation,
                           reaxpro_wrapper, emmo.hasPart)
if search_calculation:
    pretty_print(search_calculation[0])
    export_cuds(search_calculation[0], "ams_energy_landscape.ttl", format="ttl")
