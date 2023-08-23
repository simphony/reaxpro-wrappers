import os
import warnings

from setuptools import setup

PATH = os.path.dirname(__file__)
HOST = "raw.githubusercontent.com"
ROUTE = "simphony/reaxpro-framework-ontology"
VERSION: str = "v1.0.0"
REAXPRO_ONTOLOGY = "reaxpro-inferred.ttl"
URL = f"https://{HOST}/{ROUTE}/{VERSION}/{REAXPRO_ONTOLOGY}"
YML_PATH = "osp.core.ontology.docs"
REAXPRO_YML: str = f"{YML_PATH}.reaxpro"


def install_ontology():
    """Generates and installs the reaxpro ontology."""
    from tempfile import NamedTemporaryFile

    from osp.core.ontology.installation import OntologyInstallationManager

    manager = OntologyInstallationManager()

    import requests
    import yaml

    with NamedTemporaryFile(delete=False) as ontology_file:
        response = requests.get(URL)
        if response.status_code == 200:
            ontology_file.write(response.content)
            ontology_file.flush()
        else:
            message = f"""Ontology file cannot be fetched from `{URL}`.
            Status code: {response.status_code}."""
            raise RuntimeError(message)

    filepath = os.path.join(PATH, *REAXPRO_YML.split("."))

    with open(f"{filepath}.yml", "r") as file:
        content = yaml.safe_load(file)

    content["ontology_file"] = ontology_file.name

    with NamedTemporaryFile("w", suffix=".yml", delete=False) as target_yml_file:
        yaml.dump(content, target_yml_file)

    manager.install_overwrite(target_yml_file.name)

# install wrapper
setup()

# install ontology
try:
    install_ontology()
except Exception as err:
    warnings.warn(f"Ontologies were not properly installed: {err}")
