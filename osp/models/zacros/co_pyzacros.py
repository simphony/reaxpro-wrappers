from io import StringIO
import os
from pathlib import Path
import tempfile
from typing import TYPE_CHECKING, List, Union, Optional
from uuid import UUID
import warnings

from arcp import arcp_random
from pydantic import AnyUrl, BaseModel, Field, root_validator, confloat, conlist, conint
from pydantic.dataclasses import dataclass
from urllib.parse import quote, urlencode


from osp.core.cuds import Cuds
from osp.core.namespaces import emmo, crystallography
from osp.core.session import CoreSession
from osp.core.utils import export_cuds
from osp.tools.io_functions import raise_warning
from osp.core.utils import pretty_print
from osp.models.utils.general import get_upload, _get_example_json
from osp.models.multiscale.co_pt111_meso import COMolarFractionRange


STANDARD_DATA = [("9a1e6c07-840a-4182-8ed5-60212167aa4b", "lattice_input.dat")]


if TYPE_CHECKING:
    from typing import List
    from osp.core.cuds import Cuds


class Lattice(BaseModel):
    """Zacros lattice model."""

    xyz_file: Union[UUID, AnyUrl, Path] = Field(
        ...,
        description="""
        UUID of the cache-upload or url/system path to Zacros lattice input."""
    )


class ClusterExpansion(BaseModel):
    """Zacros interaction energy model."""

    name: str = Field(
       ..., description="""
       String name of the cluster.
       """
    )

    neighboring: Optional[list] = Field(
       None, description="""
       Specifies the neighboring between sites, if more than one sites appear in the graph pattern.
       It is followed by expressions structured as int1-int2 in the same line of input as the
       keyword.
       """
       )

    sites: int = Field(
       ..., description="""
       Specifies the number of sites in the graph pattern representing the cluster.
       """
    )

    lattice_state: list = Field(
       ..., description="""
       Specifies the state of each site in the graph pattern representing the cluster.
       """
    )

    site_types: Optional[list] = Field(
       None, description="""
       The site types for each of the different sites of the unit cell.
       """
    )

    graph_multiplicity: Optional[int] = Field(
        None, description="""
        The multiplicity of the pattern, namely the number of times that the exact same pattern
        will be counted for a given lattice configuration.
        """
    )

    cluster_eng: float = Field(
        ..., description="""
        The energy contribution of the pattern, given as a real number following the keyword.
        """
    )

    @root_validator
    def validate_all(cls, values):

        if isinstance(values.get("neighboring"), list):
            for neighbor_pair in values.get("neighboring"):
                if "-" not in neighbor_pair:
                    raise ValueError(
                        f'{"neighboring expressions are structured as int1-int2."}')
                else:
                    for neighbor in neighbor_pair.split("-"):
                        if not neighbor.isdigit():
                            raise ValueError(
                                f'{"neighboring must be a list of pairs of integers."}')

        sites = values.get("sites")
        if (sites == 1 and len(values.get("lattice_state")) != 3):
            raise ValueError(
                f'{"lattice_state list must be: [int1, str, int2]."}')

        if (sites > 1 and len(values.get("lattice_state")) != sites):
            raise ValueError(
                f'{"Number of lattice_state must be equal to sites."}')

        if values.get("site_types"):
            if (len(values.get("site_types")) != sites):
                raise ValueError(
                    f'{"Number of site_types must be equal to sites."}')
        return values


class Variant(BaseModel):
    """Block to define mechanism variants."""

    name: str = Field(
        ..., description="""
        String name of the variant.
        """
    )

    site_types: list = Field(
        ..., description="""
        The site types for each of the different sites of the unit cell.
        """
    )

    pre_expon: float = Field(
        ..., description="""
        Specifies the pre-exponential in the Arrhenius formula giving the rate constant of that
        elementary step. Only constant values are implemented in the BaseModel.
        """
    )

    pe_ratio: float = Field(
        ..., description="""
        This keyword gives the ratio of forward over reverse pre-exponentials and is valid only
        inside a reversible elementary step specification block.
        """
    )

    activ_eng: float = Field(
        ..., description="""
        The activation energy at the zero coverage limit.
        """
    )


class ElementaryStep(BaseModel):
    """Zacros mechanism model."""

    reversible_step: Optional[str] = Field(
        None, description="""
        This blocks defines a reversible elementary step.
        """
    )

    step: Optional[str] = Field(
        None, description="""
        This blocks defines a irreversible elementary step.
        """
    )

    gas_reacs_prods: list = Field(
        ..., description="""
        Provides information about the gas species participating in the mechanism.
        """
    )

    sites: int = Field(
        ..., description="""
        Specifies the number of sites in the graph pattern representing the elementary step being
        defined.
        """
    )

    neighboring: Optional[list] = Field(
        None, description="""
        Specifies the neighboring between sites, if more than one sites appear in the graph pattern
        representing the elementary step.
        """
    )

    initial: list = Field(
        ..., description="""
        Specifies the initial state of each site in the graph pattern.
        """
    )

    final: list = Field(
        ..., description="""
        Specifies the final state of each site in the graph pattern.
        """
    )

    variant: Optional[Variant] = Field(
        None, description="""
        Variant blocks are used to reduce repetitions in the energetics_input.dat file.
        """
    )

    pre_expon: Optional[float] = Field(
        None, description="""
        Specifies the pre-exponential in the Arrhenius formula giving the rate constant of that
        elementary step. Only constant values are implemented in the BaseModel.
        """
    )

    activ_eng: Optional[float] = Field(
        None, description="""
        The activation energy at the zero coverage limit.
        """
    )

    pe_ratio: Optional[float] = Field(
        None, description="""
        This keyword gives the ratio of forward over reverse pre-exponentials and is valid only
        inside a reversible elementary step specification block.
        """
    )

    @root_validator
    def validate_all(cls, values):

        sites = values.get("sites")
        options = [values.get("initial"),
                   values.get("final")]

        for option in options:
            if sites != len(option):
                if (sites == 1 and len(option) != 3):
                    raise ValueError(
                        f"""Length of list('{option}') should be 3, e.g. '[1, "CO*", 1]'.""")

        if values.get("variant"):
            if sites != len(values.get("variant").site_types):
                raise ValueError(
                    f"""Length of 'variant/site_types': '{values.get("variant").site_types}' should
                    be equal to 'sites', '{sites}'""")
        else:

            if values.get("pre_expon") is None:
                raise ValueError(
                    f'{"""Missing pre_expon in elementary step."""}')

            if values.get("activ_eng") is None:
                raise ValueError(
                    f'{"""Missing activ_eng in elementary step."""}')

        return values


class Temperature(BaseModel):
    """Temperature base model."""

    value: confloat(gt=0, allow_inf_nan=False) = Field(
        ..., description="""
        Temperature (K) under which the system is simulated.
        """
    )

    ramp: Optional[list] = Field(
        None, description="""
        'ramp': [real1, real2] specifies a temperature ramp
        where real1 is the initial temperature (K) and real2 is the rate of change (K/s).
        If real2 is positive, temperature programmed desorption can be simulated.
        Negative values of real2 can be used for simulated annealing calculations."""
    )

    @root_validator
    def validate_all(cls, values):
        if values.get("ramp"):
            if len(values.get("ramp")) != 2:
                raise ValueError(
                    """List of 'temperature' ramp values must have two fields,
                    [temp, rate_of_change]."""
                )
        return values


@dataclass
class SimulationSettings:
    """Initiates simulation settings from dict."""

    random_seed: int = Field(
        ..., description="""
        The integer seed of the random number generator.
        """
    )

    temperature: Temperature = Field(
        ..., description="""
        Base model temperature (K) under which the system is simulated.
        """
    )

    pressure: confloat(gt=0, allow_inf_nan=False) = Field(
        ..., description="""
        The pressure (bar) under which the system is simulated.
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

    gas_energies: conlist(item_type=confloat(allow_inf_nan=False)) = Field(
        ..., description="""
        A list[str1, str2, ...] with the total energies (eV) of the gas species.
        There should be as many reals following this keyword as the number of gas species
        specified with keyword n_gas_species. The ordering of these values should be consistent
        with the order used in keyword gas_specs_names.
        """
    )

    gas_molec_weights: conlist(item_type=confloat(ge=0, allow_inf_nan=False)) = Field(
        ..., description="""
        A list[str1, str2, ...] with the molecular weights (amu) of the gas species.
        There should be as many reals following the keyword as the number of gas species specified
        with keyword n_gas_species. Note: at present these values are not used in the code. This
        feature is there for future development.
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

    n_surf_species: int = Field(
        ..., description="""
        The number of surface species in the chemistry.
        """
    )

    surf_specs_names: conlist(item_type=str) = Field(
       ..., description="""
        surf_specs_names is a list[str1, str2, ...] with the names of the surface species.
        There should be as many strings following the keyword as the number of surface species
        specified with keyword n_surf_species."
        """
    )

    surf_specs_dent: conlist(item_type=int) = Field(
       ..., description="""
        A list[int1, int2, ...] with the number of dentates of the surface species,
        specifying the number of sites each species binds to.
        """
    )

    snapshots: list = Field(
       ..., description="""
        Determines how often snapshots of the lattice state will be written to output file
        history_output.txt.
        """
    )

    process_statistics: list = Field(
       ..., description="""
        Determines how often statistical information about the occurrence of elementary
        events will be written to output file procstat_output.txt
        """
    )

    species_numbers: list = Field(
       ..., description="""
        Determines how often information about the number of gas and surface species,
        as well as the energy of the current lattice configuration) will be written
        to specnum_output.txt
        """
    )

    event_report: str = Field(
       ..., description="""
       Controls event reporting behavior. Expression expr can be 'on' or 'off'.
       """
    )

    max_steps: confloat(ge=0, allow_inf_nan=True) = Field(
       ..., description="""
       The maximum number of KMC steps to be simulated. This keyword defines a stopping criterion.
       """
    )

    max_time: confloat(ge=0, allow_inf_nan=True) = Field(
       ..., description="""
       The maximum allowed simulated time interval (time ranges from 0.0 to the maximum time in a
       simulation).
       """
    )

    wall_time: conint(gt=0) = Field(
       ..., description="""
       The maximum allowed real-world time in seconds that the simulation can be left running.
       The code has an internal "stopwatch” and will exit normally once the wall time has been
       exceeded.
       """
    )

    @root_validator
    def validate_all(cls, values):

        gas_variables = ["gas_specs_names",
                         "gas_energies",
                         "gas_molec_weights",
                         "gas_molar_fracs"]

        surf_variables = ["surf_specs_names",
                          "surf_specs_dent"]

        recoding_quantities = ["snapshots",
                               "process_statistics",
                               "species_numbers"]

        recording_options = ["off",
                             "on time",
                             "on event",
                             "on elemevent",
                             "on logtime",
                             "on realtime"]

        # Validate gas species:
        for variable in gas_variables:
            if len(values.get(variable)) != values.get("n_gas_species"):
                raise ValueError(
                    f"""Length of '{variable}' list must be equal to 'n_gas_species'.""")

        # Validate adsorbed species:
        for variable in surf_variables:
            if len(values.get(variable)) != values.get("n_surf_species"):
                raise ValueError(
                    f"""Length of '{variable}' list must be equal to 'n_surf_species'.""")

        # Validate recording options:
        for quantity in recoding_quantities:
            if values[quantity][0] not in recording_options:
                raise ValueError(
                    f"""Wrong recording option in '{quantity}', use one of the following:
                                "off", "on time", "on event", "on elemevent",
                                "on logtime", "on realtime" """)

            if values[quantity][0] == "on event":
                if len(values[quantity]) == 1:
                    values[quantity].append(1)

            if len(values[quantity]) != 2:
                raise ValueError(
                    f"""Wrong length of '{quantity}' list, you should add the
                    '{values[quantity][0]}' frequency
                    recording.""")

            elif values[quantity][0] == "on logtime":
                if len(values[quantity]) != 3:
                    raise ValueError(
                        f"""Wrong length of '{quantity}' list, you should add the
                        '{values[quantity][0]}' the initial
                        time and logarithmic constant.""")

        # Validate event_report:
        recording_options = ["on", "off"]

        if values.get("event_report") not in recording_options:
            raise ValueError(f'{"""Wrong event_report option, choose on or off"""}')

        return values


@dataclass
class COpyZacrosModel:
    """Initiates a mesoscopic calculation from dict."""

    simulation_input: SimulationSettings = Field(
        ..., description="Zacros simulation setting model.")

    lattice_input: Lattice = Field(
        ..., description="Zacros lattice model.")

    energetics_input: List[ClusterExpansion] = Field(
        ..., description="Zacros interaction energy model.")

    mechanism_input: List[ElementaryStep] = Field(
        ..., description="List of chemical reactions.")

    adp: Optional[COMolarFractionRange] = Field(
        None,
        description="""Molar fractions of CO
        for the adaptive design procedure"""
        )

    @root_validator
    def validate_all(cls, values):

        adsorbed_species_from_settings = values.get("simulation_input").surf_specs_names
        adsorbed_species_from_clusters = []

        clusters = values.get("energetics_input")
        for cluster in clusters:
            for lattice in cluster.lattice_state:
                if isinstance(lattice, list):
                    adsorbed_species_from_clusters.extend(species for species in lattice)
                else:
                    adsorbed_species_from_clusters.append(lattice)
        adsorbed_species_from_clusters = list(dict.fromkeys([
            species for species in adsorbed_species_from_clusters if isinstance(species, str)]))

        for species in adsorbed_species_from_settings:
            if species not in adsorbed_species_from_clusters:
                raise ValueError(
                                f"""'{species}' not found in 'lattice_state'
                                of 'energetics_input'.""")

        for species in adsorbed_species_from_clusters:
            if species not in adsorbed_species_from_settings:
                raise ValueError(
                                f"""'{species}' not found in 'surf_specs_names'
                                of 'simulation_input'.""")

        # Validation of gas_reacs_prods from mechanism vs gas_specs_names from settings:

        gas_species_from_settings = values.get("simulation_input").gas_specs_names
        gas_species_from_mechanism = []

        steps = values.get("mechanism_input")

        for step in steps:
            for species in step.gas_reacs_prods:
                # more than one gas product:
                if isinstance(species, list):
                    for item in species:
                        if isinstance(item, str):
                            gas_species_from_mechanism.append(item)
                elif isinstance(species, str):
                    gas_species_from_mechanism.append(species)

        for species in gas_species_from_settings:
            if species not in gas_species_from_mechanism:
                raise ValueError(
                      f"""'{species}' not found in 'gas_reacs_prods' of 'mechanism_input'.""")

        for species in gas_species_from_mechanism:
            if species not in gas_species_from_settings:
                raise ValueError(
                      f"""'{species}' not found in 'gas_specs_names' of 'simulation_input'.""")

        return values

    def __post_init_post_parse__(self):
        with CoreSession() as session:

            if not self.adp:

                calculation = emmo.MesoscopicCalculation()
                calculation.add(*self._make_model(), rel=emmo.hasInput)

            else:
                workflow = emmo.Workflow()
                calculation = emmo.MesoscopicCalculation()
                calculation.add(*self._make_model(), rel=emmo.hasInput)
                apd = emmo.AdaptiveDesignProcedure()
                apd.add(self.adp.cuds, rel=emmo.hasInput)
                workflow.add(calculation, rel=emmo.hasSpatialFirst)
                calculation.add(apd, rel=emmo.hasSpatialNext)
                workflow.add(apd, rel=emmo.hasSpatialLast)

        file = tempfile.NamedTemporaryFile(suffix=".ttl", delete=False)
        export_cuds(session, file.name)
        self._file = file.name
        try:
            self._uuid = get_upload(file)
        except Exception as error:
            self._uuid = None
            message = message = f"The graph of the model could not be stored at the minio-instance: {error.args}"
            warnings.warn(message)
        self._session = session

    def _make_arcp(self, *args, **kwargs):
        if kwargs.get("query"):
            query = kwargs.pop("query")
            for key, value in query.items():
                query[key] = [quote(".".join(item)) for item in value]
            query = urlencode(query, doseq=True)
        return arcp_random(*args, **kwargs, query=query)


    def _make_model(self) -> "List[Cuds]":
        """Construction of model CUDS graph."""
        cuds = []

        # simulation_input.dat
        cuds.extend(self._make_settings())

        # lattice_input.dat
        cuds.append(self._make_lattice())

        # energetics_input.dat
        cuds.extend(self._make_energetics())

        # mechanism_input.dat
        cuds.append(self._make_mechanism())

        return cuds

    def _make_settings(self) -> "List[Cuds]":
        """Construction of setting CUDS graph."""

        cuds_list = []
        cuds_list.append(self._make_random_seed())
        cuds_list.append(self._make_temperature())
        cuds_list.append(self._make_pressure())
        cuds_list.extend(self._make_gas_species())
        cuds_list.extend(self._make_surf_species())
        cuds_list.extend(self._make_recordings())
        cuds_list.append(self._make_event_report())
        cuds_list.append(self._make_max_steps())
        cuds_list.append(self._make_max_time())
        cuds_list.append(self._make_wall_time())

        return cuds_list

    def _make_random_seed(self) -> "Cuds":
        """Add random_seed to CUDS graph."""

        iri = self._make_arcp("random_seed", query=dict(jsonpath=[["random_seed"]]))
        random_seed = emmo.RandomSeed(hasNumericalData=self.simulation_input.random_seed, iri=iri)

        return random_seed

    def _make_temperature(self) -> "Cuds":
        """Add temperature to CUDS graph."""

        if self.simulation_input.temperature.ramp:

            raise_warning(file=os.path.basename(__file__),
                          function=self._make_temperature.__name__,
                          message="Ramp temperature excluded from knowledge graph!")

        iri = self._make_arcp("temperature", query=dict(jsonpath=[["temperature", "value"]]))
        temperature_float = emmo.Real(
            hasNumericalData=self.simulation_input.temperature.value, iri=iri)
        temperature_unit = emmo.Kelvin(hasSymbolData='K')
        temperature = emmo.ThermodynamicTemperature()
        temperature.add(temperature_unit, rel=emmo.hasReferenceUnit)
        temperature.add(temperature_float, rel=emmo.hasQuantityValue)

        return temperature

    def _make_pressure(self) -> "Cuds":
        """Add pressure to CUDS graph."""

        iri = self._make_arcp("pressure", query=dict(jsonpath=[["pressure"]]))
        pressure_float = emmo.Real(
            hasNumericalData=self.simulation_input.pressure * 101325, iri=iri)
        # Bar unit not implemented in emmo.

        pressure_unit = emmo.Pascal(hasSymbolData='Pa')
        pressure = emmo.Pressure()
        pressure.add(pressure_unit, rel=emmo.hasReferenceUnit)
        pressure.add(pressure_float, rel=emmo.hasQuantityValue)

        return pressure

    def _make_gas_species(self) -> "List[Cuds]":
        """Add gas species to CUDS graph"""

        cuds = []
        for species in range(self.simulation_input.n_gas_species):

            gas_species = emmo.GasSpecies()
            formation_energy = emmo.FormationEnergy()
            molar_fraction = emmo.AmountFraction()
            gas_molec_weight = emmo.MolecularWeight()

            iri = self._make_arcp("gas_specs_names",
                            query=dict(jsonpath=[["gas_specs_names", str(species)]]))
            species_name = emmo.ChemicalElement(
                hasSymbolData=self.simulation_input.gas_specs_names[species], iri=iri)

            iri = self._make_arcp("gas_energies", query=dict(jsonpath=[["gas_energies", str(species)]]))
            formation_energy_value = emmo.Real(
                hasNumericalData=self.simulation_input.gas_energies[species], iri=iri)
            formation_energy.add(formation_energy_value, rel=emmo.hasQuantityValue)

            iri = self._make_arcp(
                "gas_molec_weights", query=dict(jsonpath=[["gas_molec_weights", str(species)]]))
            gas_molec_weight_float = emmo.Real(
                hasNumericalData=self.simulation_input.gas_molec_weights[species], iri=iri)
            gas_molec_weight.add(gas_molec_weight_float, rel=emmo.hasQuantityValue)

            iri = self._make_arcp("gas_molar_fracs",
                            query=dict(jsonpath=[["gas_molar_fracs", str(species)]]))
            molar_fraction_float = emmo.Real(
                hasNumericalData=self.simulation_input.gas_molar_fracs[species], iri=iri)
            molar_fraction.add(molar_fraction_float, rel=emmo.hasQuantityValue)

            gas_species.add(species_name,
                            formation_energy,
                            gas_molec_weight, molar_fraction, rel=emmo.hasProperty)
            cuds.append(gas_species)

        return cuds

    def _make_surf_species(self) -> "List[Cuds]":
        """Add surface-adsorbed species to CUDS graph"""

        cuds_list = []

        for species in range(self.simulation_input.n_surf_species):

            surf_species = emmo.AdsorbedSpecies()
            denticity = emmo.Denticity()
            iri = self._make_arcp("surf_specs_names",
                            query=dict(jsonpath=[["surf_specs_names", str(species)]]))
            species_name = emmo.ChemicalElement(
                hasSymbolData=self.simulation_input.surf_specs_names[species], iri=iri)
            iri = self._make_arcp("surf_specs_dent",
                            query=dict(jsonpath=[["surf_specs_dent", str(species)]]))
            denticity_number = emmo.Real(
                hasNumericalData=self.simulation_input.surf_specs_dent[species], iri=iri)
            denticity.add(denticity_number, rel=emmo.hasQuantityValue)
            surf_species.add(species_name, denticity, rel=emmo.hasSpatialDirectPart)
            cuds_list.append(surf_species)

        return cuds_list

    def _make_recordings(self) -> "List[Cuds]":
        """Add output recording options to CUDS graph"""

        cuds_list = []

        recoding_quantities = ["snapshots",
                               "process_statistics",
                               "species_numbers"]

        for quantity in recoding_quantities:
            cuds_quantity = self._evaluate_recording_quantity(quantity)
            cuds_list.append(cuds_quantity)

        return cuds_list

    def _evaluate_recording_quantity(self, quantity: str) -> Cuds:
        """Evaluate recording option for a given output type"""

        if quantity == "snapshots":
            recording_option = self.simulation_input.snapshots[0]
            iri = self._make_arcp("snapshots", query=dict(jsonpath=[["snapshot", str(0)]]))
            cuds_object = emmo.Snapshots(hasSymbolData=recording_option, iri=iri)
            if recording_option != "off":
                cuds_object.add(
                    self._evaluate_recording_combination(self.simulation_input.snapshots,
                                                         "snapshots"), rel=emmo.hasSpatialPart)

        if quantity == "process_statistics":
            recording_option = self.simulation_input.process_statistics[0]
            iri = self._make_arcp("process_statistics",
                            query=dict(jsonpath=[["process_statistics", str(0)]]))
            cuds_object = emmo.ProcessStatistics(hasSymbolData=recording_option, iri=iri)
            if recording_option != "off":
                cuds_object.add(
                    self._evaluate_recording_combination(
                        self.simulation_input.process_statistics,
                        "process_statistics"), rel=emmo.hasSpatialPart)

        if quantity == "species_numbers":
            recording_option = self.simulation_input.species_numbers[0]
            iri = self._make_arcp("species_numbers", query=dict(jsonpath=[["species_numbers", str(0)]]))
            cuds_object = emmo.SpeciesNumbers(hasSymbolData=recording_option)
            if recording_option != "off":
                cuds_object.add(
                    self._evaluate_recording_combination(
                        self.simulation_input.species_numbers,
                        "species_numbers"), rel=emmo.hasSpatialPart)

        return cuds_object

    def _evaluate_recording_combination(self, recording_options: list, keyword: str) -> Cuds:
        """Evaluate recording freq. and type given a list of recording options"""

        option = recording_options[0]

        if option == "on logtime":
            cuds_object = emmo.Array()
            iri = self._make_arcp(keyword, query=dict(jsonpath=[[keyword, str(1)]]))
            time_float_1 = emmo.Real(hasNumericalData=recording_options[1], iri=iri)
            iri = self._make_arcp(keyword, query=dict(jsonpath=[[keyword, str(2)]]))
            time_float_2 = emmo.Real(hasNumericalData=recording_options[2], iri=iri)
            cuds_object.add(time_float_1, time_float_2, rel=emmo.hasSpatialPart)

        else:
            iri = self._make_arcp(keyword, query=dict(jsonpath=[[keyword, str(3)]]))
            cuds_object = emmo.Real(hasNumericalData=recording_options[1], iri=iri)

        return cuds_object

    def _make_event_report(self) -> Cuds:
        """Creates CUDS graph of the event on/off switch"""

        iri = self._make_arcp("event_report", query=dict(jsonpath=[["event_report"]]))
        if self.simulation_input.event_report == "off":
            cuds_object = emmo.EventReport(hasNumericalData=False, iri=iri)

        elif self.simulation_input.event_report == "on":
            cuds_object = emmo.EventReport(hasNumericalData=False, iri=iri)

        return cuds_object

    def _make_max_steps(self) -> Cuds:
        """Creates CUDS graph for the max steps option"""

        iri = self._make_arcp("max_steps", query=dict(jsonpath=[["max_steps"]]))
        cuds_object = emmo.MaximumSteps()
        step_float = emmo.Real(hasNumericalData=self.simulation_input.max_steps, iri=iri)
        cuds_object.add(step_float, rel=emmo.hasPart)

        return cuds_object

    def _make_max_time(self) -> Cuds:
        """Creates CUDS graph for the max time option"""

        iri = self._make_arcp("max_time", query=dict(jsonpath=[["max_time"]]))
        cuds_object = emmo.MaximumTime()
        step_float = emmo.Real(hasNumericalData=self.simulation_input.max_time, iri=iri)
        cuds_object.add(step_float, rel=emmo.hasPart)

        return cuds_object

    def _make_wall_time(self) -> Cuds:
        """Creates CUDS graph for the wall time option"""

        iri = self._make_arcp("wall_time", query=dict(jsonpath=[["wall_time"]]))
        cuds_object = emmo.WallTime()
        step_float = emmo.Real(hasNumericalData=self.simulation_input.wall_time, iri=iri)
        cuds_object.add(step_float, rel=emmo.hasPart)

        return cuds_object

    def _make_lattice(self) -> Cuds:
        """Creates CUDS graph for the lattice input"""
        xyz_file = self.lattice_input.xyz_file
        if isinstance(xyz_file, UUID):
            return crystallography.UnitCell(uid=xyz_file)
        elif isinstance(xyz_file, AnyUrl):
            return crystallography.UnitCell(iri=str(xyz_file))
        elif isinstance(xyz_file, Path):
            if not "file:" in str(xyz_file):
                iri = f"file://{xyz_file.as_posix()}"
            else:
                iri = xyz_file.as_posix()
            return crystallography.UnitCell(iri=iri)

    def _make_energetics(self) -> "List[Cuds]":
        """Creates CUDS graph for the energy cluster input"""

        cuds_list = []

        for cluster in self.energetics_input:

            num_of_sites = cluster.sites

            iri = self._make_arcp(
                "energetics_input",
                query=dict(jsonpath=[["energetics_input", "name", cluster.name]]))
            cuds_object = emmo.ClusterExpansion(iri=iri)

            # neighboring
            if isinstance(cluster.neighboring, list):

                for neighbor_pair in cluster.neighboring:

                    neighbor_site_1 = emmo.Real(hasNumericalData=neighbor_pair.split("-")[0])
                    neighbor_site_2 = emmo.Real(hasNumericalData=neighbor_pair.split("-")[1])
                    iri = self._make_arcp("energetics_input",
                                    query=dict(jsonpath=[["energetics_input", "neighboring"]]))
                    neighbors = emmo.Neighboring(iri=iri)
                    neighbors.add(neighbor_site_1, rel=emmo.hasSpatialFirst)
                    neighbors.add(neighbor_site_2, rel=emmo.hasSpatialNext)

                cuds_object.add(neighbors, rel=emmo.hasSpatialPart)

            # lattice_state
            if num_of_sites == 1:

                iri = self._make_arcp("energetics_input",
                                query=dict(jsonpath=[["energetics_input", "lattice_state"]]))
                lattice_state = emmo.LatticeState(iri=iri)
                lattice_state.add(self._get_lattice_state_from_list(cluster.lattice_state),
                                  rel=emmo.hasSpatialDirectPart)
                cuds_object.add(lattice_state, rel=emmo.hasPart)

            else:

                for index, state in enumerate(range(num_of_sites)):

                    iri = self._make_arcp("energetics_input",
                                    query=dict(jsonpath=[[
                                        "energetics_input", "lattice_state", str(index)]]))
                    lattice_state = emmo.LatticeState(iri=iri)
                    lattice_state.add(
                        self._get_lattice_state_from_list(
                            cluster.lattice_state[state], index=index),
                        rel=emmo.hasSpatialDirectPart)
                    cuds_object.add(lattice_state, rel=emmo.hasPart)

            # site_types
            if isinstance(cluster.site_types, list):

                for index, site_type in enumerate(range(num_of_sites)):

                    iri = self._make_arcp("energetics_input",
                                    query=dict(jsonpath=[["energetics_input",
                                                          "site_types", str(index)]]))
                    site_type = crystallography.Site()
                    atom_symbol = emmo.ChemicalElement(hasSymbolData=cluster.site_types[index], iri=iri)
                    site_type.add(atom_symbol, rel=emmo.hasSpatialDirectPart)
                    cuds_object.add(site_type, rel=emmo.hasSpatialPart)

            # cluster_eng
            iri = self._make_arcp("energetics_input",
                            query=dict(jsonpath=[["energetics_input", "cluster_eng"]]))
            cluster_energy = emmo.ClusterEnergy(iri=iri)
            cluster_energy_value = emmo.Real(hasNumericalData=cluster.cluster_eng)
            cluster_energy.add(cluster_energy_value, rel=emmo.hasQuantityValue)
            cuds_object.add(cluster_energy, rel=emmo.hasSpatialPart)

            cuds_list.append(cuds_object)

        return cuds_list

    def _get_lattice_state_from_list(self, lattice_state: list, index=0) -> Cuds:
        """Generates graph for a given list of a lattice_state."""

        cuds_object = emmo.AdsorbedSpecies()
        if index != 0:
            index = index + 1

        iri = self._make_arcp("energetics_input",
                        query=dict(jsonpath=[["energetics_input", "lattice_state", str(index)]]))
        species_number = emmo.Real(hasNumericalData=lattice_state[0], iri=iri)

        iri = self._make_arcp("energetics_input",
                        query=dict(jsonpath=[["energetics_input", "lattice_state", str(index+1)]]))
        atom_symbol = emmo.ChemicalElement(hasSymbolData=lattice_state[1], iri=iri)

        iri = self._make_arcp("energetics_input",
                        query=dict(jsonpath=[["energetics_input", "lattice_state", str(index+2)]]))
        denticity_number = emmo.Real(hasNumericalData=lattice_state[2], iri=iri)
        denticity = emmo.Denticity()
        denticity.add(denticity_number, rel=emmo.hasQuantityValue)

        cuds_object.add(species_number, atom_symbol, denticity, rel=emmo.hasSpatialDirectPart)
        return cuds_object

    def _make_mechanism(self) -> Cuds:
        """Construction of mechanism CUDS graph."""

        cuds_object = emmo.ChemicalReactionMechanism()

        for step in self.mechanism_input:

            if step.reversible_step:

                iri = self._make_arcp("mechanism_input", query=dict(
                    jsonpath=[["mechanism_input", "reversible_step", step.reversible_step]]))
                reaction = emmo.ReversibleReactionEquation(iri=iri)

            elif step.step:

                iri = self._make_arcp("mechanism_input", query=dict(
                    jsonpath=[["mechanism_input", "step", step.step]]))
                reaction = emmo.IrreversibleReactionEquation(iri=iri)

            # gas_reacs_prods
            for index, species in enumerate(step.gas_reacs_prods):

                # more than one gas product:
                if isinstance(species, list):
                    reaction.add(self._get_gas_species_from_list(
                        species, index=index), rel=emmo.hasSpatialPart)
                # just one:
                elif isinstance(species, str):
                    reaction.add(self._get_gas_species_from_list(step.gas_reacs_prods),
                                 rel=emmo.hasSpatialPart)

            # neighboring
            if step.neighboring:
                if len(step.neighboring) > 1:
                    for index, neighbor in enumerate(step.neighboring):
                        reaction.add(self._get_neighbor_from_string(neighbor, index=index),
                                     rel=emmo.hasSpatialPart)
                else:
                    reaction.add(self._get_neighbor_from_string(step.neighboring[0]),
                                 rel=emmo.hasSpatialPart)

            # initial
            if step.sites == 1:
                reactant = emmo.ChemicalReactionEquationReactant()
                iri = self._make_arcp("mechanism_input",
                                query=dict(jsonpath=[["mechanism_input", "initial", str(1)]]))
                atom_symbol = emmo.ChemicalElement(hasSymbolData=step.initial[1], iri=iri)
                iri = self._make_arcp("mechanism_input",
                                query=dict(jsonpath=[["mechanism_input", "initial", str(2)]]))
                coefficient = emmo.StoichiometricCoefficient(iri=iri)
                coefficient_value = emmo.Real(hasNumericalData=step.initial[2])
                coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                reactant.add(
                    atom_symbol, coefficient, rel=emmo.hasSpatialDirectPart
                )
                reaction.add(reactant, rel=emmo.hasSpatialPart)
            else:

                for site in range(step.sites):

                    reactant = emmo.ChemicalReactionEquationReactant()
                    iri = self._make_arcp("mechanism_input",
                                    query=dict(
                                        jsonpath=[[
                                            "mechanism_input", "initial", str(site), str(1)]]))
                    atom_symbol = emmo.ChemicalElement(
                        hasSymbolData=(step.initial[site])[1], iri=iri)
                    iri = self._make_arcp("mechanism_input",
                                    query=dict(
                                        jsonpath=[[
                                            "mechanism_input", "initial", str(site), str(2)]]))
                    coefficient = emmo.StoichiometricCoefficient(iri=iri)
                    coefficient_value = emmo.Real(hasNumericalData=(step.initial[site])[2])
                    coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                    reactant.add(
                        atom_symbol, coefficient, rel=emmo.hasSpatialDirectPart
                    )
                    reaction.add(reactant, rel=emmo.hasSpatialPart)

            # final
            if step.sites == 1:

                reactant = emmo.ChemicalReactionEquationProduct()
                iri = self._make_arcp("mechanism_input",
                                query=dict(jsonpath=[["mechanism_input", "final", str(1)]]))
                atom_symbol = emmo.ChemicalElement(hasSymbolData=step.final[1], iri=iri)
                iri = self._make_arcp("mechanism_input",
                                query=dict(jsonpath=[["mechanism_input", "final", str(2)]]))
                coefficient = emmo.StoichiometricCoefficient(iri=iri)
                coefficient_value = emmo.Real(hasNumericalData=step.final[2])
                coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                reactant.add(
                    atom_symbol, coefficient, rel=emmo.hasSpatialDirectPart
                )
                reaction.add(reactant, rel=emmo.hasSpatialPart)

            else:

                for site in range(step.sites):

                    reactant = emmo.ChemicalReactionEquationProduct()
                    iri = self._make_arcp(
                        "mechanism_input",
                        query=dict(jsonpath=[["mechanism_input", "final", str(site), str(1)]]))
                    atom_symbol = emmo.ChemicalElement(hasSymbolData=(step.final[site])[1], iri=iri)
                    iri = self._make_arcp("mechanism_input",
                                    query=dict(jsonpath=[["mechanism_input", "final",
                                                          str(site), str(2)]]))
                    coefficient = emmo.StoichiometricCoefficient(iri=iri)
                    coefficient_value = emmo.Real(hasNumericalData=(step.final[site])[2])
                    coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                    reactant.add(
                        atom_symbol, coefficient, rel=emmo.hasSpatialDirectPart
                    )
                    reaction.add(reactant, rel=emmo.hasSpatialPart)

            # site_types
            if step.variant:

                for site in range(step.sites):

                    iri = self._make_arcp(
                        "mechanism_input",
                        query=dict(jsonpath=[["mechanism_input", "variant", "site_types",
                                              str(step.variant.site_types[site])]]))
                    site_type = crystallography.Site(iri=iri)
                    atom_symbol = emmo.ChemicalElement(hasSymbolData=step.variant.site_types[site])
                    site_type.add(atom_symbol, rel=emmo.hasSpatialDirectPart)
                    reaction.add(site_type, rel=emmo.hasSpatialPart)

                # pre_expon
                iri = self._make_arcp(
                    "mechanism_input", query=dict(jsonpath=[["mechanism_input", "variant",
                                                             "pre_expon"]]))
                coefficient = emmo.ArrheniusCoefficient(iri=iri)
                coefficient_value = emmo.Real(hasNumericalData=step.variant.pre_expon)
                coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                reaction.add(coefficient, rel=emmo.hasSpatialPart)

                # pe_ratio
                iri = self._make_arcp(
                    "mechanism_input", query=dict(jsonpath=[["mechanism_input", "variant",
                                                             "pe_ratio"]]))
                coefficient = emmo.Constant(iri=iri)
                coefficient_value = emmo.Real(hasNumericalData=step.variant.pe_ratio)
                coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                reaction.add(coefficient, rel=emmo.hasSpatialPart)

                # activ_eng
                iri = self._make_arcp("mechanism_input",
                                query=dict(jsonpath=[["mechanism_input", "variant", "activ_eng"]]))
                coefficient = emmo.ActivationEnergy(iri=iri)
                coefficient_value = emmo.Real(hasNumericalData=step.variant.activ_eng)
                coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                reaction.add(coefficient, rel=emmo.hasSpatialPart)

                cuds_object.add(reaction, rel=emmo.hasPart)

            else:
                # if not Variant
                iri = self._make_arcp(
                    "mechanism_input", query=dict(jsonpath=[["mechanism_input", "pre_expon"]]))
                coefficient = emmo.ArrheniusCoefficient(iri=iri)
                coefficient_value = emmo.Real(hasNumericalData=step.pre_expon)
                coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                reaction.add(coefficient, rel=emmo.hasSpatialPart)

                # pe_ratio (optional)
                if step.pe_ratio:
                    iri = self._make_arcp(
                        "mechanism_input", query=dict(jsonpath=[["mechanism_input",
                                                                 "pe_ratio"]]))
                    coefficient = emmo.Constant(iri=iri)
                    coefficient_value = emmo.Real(hasNumericalData=step.pe_ratio)
                    coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                    reaction.add(coefficient, rel=emmo.hasSpatialPart)

                # activ_eng
                iri = self._make_arcp("mechanism_input",
                                query=dict(jsonpath=[["mechanism_input", "activ_eng"]]))
                coefficient = emmo.ActivationEnergy(iri=iri)
                coefficient_value = emmo.Real(hasNumericalData=step.activ_eng)
                coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)
                reaction.add(coefficient, rel=emmo.hasSpatialPart)
                cuds_object.add(reaction, rel=emmo.hasPart)

        return cuds_object

    def _get_gas_species_from_list(self, species_pair: list, index=0) -> Cuds:
        """Generate a (reactant/product) gas species from a list("species", "float")"""

        if index != 0:
            index = index + 1

        if float(species_pair[1]) < 0:
            cuds_object = emmo.GasReactantSpecies()

        elif float(species_pair[1]) > 0:
            cuds_object = emmo.GasProductSpecies()

        iri = self._make_arcp("simulation_input", query=dict(
            jsonpath=[["simulation_input", "gas_energies",
                      str(self._get_formation_energy_from_label(species_pair[0])[1])]]))
        formation_energy = emmo.FormationEnergy(iri=iri)
        formation_energy_value = emmo.Real(
            hasNumericalData=self._get_formation_energy_from_label(species_pair[0])[0])
        formation_energy.add(formation_energy_value, rel=emmo.hasQuantityValue)

        iri = self._make_arcp("mechanism_input", query=dict(
            jsonpath=[["mechanism_input", "gas_reacs_prods", str(index)]]))
        atom_symbol = emmo.ChemicalElement(hasSymbolData=species_pair[0], iri=iri)

        iri = self._make_arcp("mechanism_input", query=dict(
            jsonpath=[["mechanism_input", "gas_reacs_prods", str(index+1)]]))
        coefficient = emmo.StoichiometricCoefficient(iri=iri)
        coefficient_value = emmo.Real(hasNumericalData=species_pair[1])
        coefficient.add(coefficient_value, rel=emmo.hasQuantityValue)

        cuds_object.add(atom_symbol, coefficient, formation_energy,
                        rel=emmo.hasSpatialDirectPart)

        return cuds_object

    def _get_formation_energy_from_label(self, species_label: str) -> list:
        """Given a label of a gas species, search its formation energy in 'simulation_input'"""

        index = self.simulation_input.gas_specs_names.index(species_label)
        formation_energy = self.simulation_input.gas_energies[index]

        return [formation_energy, index]

    def _get_neighbor_from_string(self, neighbor: str, index=0) -> Cuds:
        """Generates a neighbor graph from a string e.g. [1-2]"""

        num_neighbor_sites = neighbor.split("-")
        cuds_object = emmo.Neighboring()

        for index2, neighbor in enumerate(num_neighbor_sites):

            iri = self._make_arcp("mechanism_input",
                            query=dict(
                                jsonpath=[[
                                    "mechanism_input", "neighboring", str(index), str(index2)]]))

            neighbor_site = emmo.Real(hasNumericalData=neighbor, iri=iri)
            if index2 == 0:
                cuds_object.add(neighbor_site, rel=emmo.hasSpatialFirst)
            else:
                cuds_object.add(neighbor_site, rel=emmo.hasSpatialNext)

        return cuds_object

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
        """Pydantic Config for COpyZacrosModel"""

        schema_extra = {
            "example": _get_example_json("Ziff-Gulari-Barshad-model-variants.json", STANDARD_DATA)
        }
