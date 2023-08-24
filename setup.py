"""Setup class for ReaxPro-wrappers"""
from setuptools import setup
import warnings

setup()

VERSION: str = "v1.0.0"
BASE = "https://raw.githubusercontent.com/simphony/reaxpro-framework-ontology"
URL = f"{BASE}/{VERSION}/reaxpro-inferred.ttl"
CONTENT = f"""identifier: reaxpro
ontology_file: '{URL}'
reference_by_label: True
format: ttl
namespaces:
    emmo: http://emmo.info/emmo
    crystallography: http://emmo.info/domain-crystallography/crystallography
active_relationships:
    - http://emmo.info/emmo#EMMO_8c898653_1118_4682_9bbf_6cc334d16a99
    - http://emmo.info/emmo#EMMO_60577dea_9019_4537_ac41_80b0fb563d41
    - http://emmo.info/emmo#EMMO_8785be5a_2493_4b12_8f39_31907ab11748
    - http://emmo.info/emmo#EMMO_52d08d7d_e9e4_43e5_8508_d353e7e3a23a
default_relationship: http://emmo.info/emmo#EMMO_17e27c22_37e1_468c_9dd7_95e137f73e7f
"""


def install_ontology():
    """Generates and installs the reaxpro ontology."""
    from tempfile import NamedTemporaryFile

    from osp.core.ontology.installation import OntologyInstallationManager

    with NamedTemporaryFile("w", suffix=".yml", delete=False) as yml_file:
        yml_file.write(CONTENT)

    manager = OntologyInstallationManager()
    manager.install_overwrite(yml_file.name)


# install ontology
try:
    install_ontology()
except Exception as err:
    warnings.warn(f"Ontologies were not properly installed: {err}")
