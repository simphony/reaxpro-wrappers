"""Conftest for simphony-catalytic"""
import os
import tempfile
from typing import TextIO, Union
from uuid import uuid4

import pytest
from minio.error import MinioException

PATH = os.path.dirname(__file__)
TEST_PATH = os.path.join(PATH, "test_files")


def return_file(key: str, as_file=False) -> str:  # pylint: disable=unused-argument
    """Mocker function for returning test files."""
    tempdir = tempfile.gettempdir()
    path = os.path.join(tempdir, f"{key}.xyz")
    if not os.path.exists(path):
        raise MinioException(f"Key `{key}` does not exist in cache.")
    return path


def store_file(file: Union[TextIO, str], uuid: str = None) -> str:
    """Mocker for uploading files."""
    tempdir = tempfile.gettempdir()
    if not uuid:
        uuid = str(uuid4())
    path = os.path.join(tempdir, f"{uuid}.xyz")
    with open(path, "wb") as temp:
        temp.write(file.read())
    return uuid


@pytest.fixture
def prepare_env():
    os.environ["CELERY_MINIO_USER"] = "foo"
    os.environ["CELERY_MINIO_PASSWORD"] = "bar"


@pytest.fixture(autouse=True)
@pytest.mark.usefixtures("mocker")
def mock_download(mocker, prepare_env):
    """Sending request to mock-minio."""
    function = "osp.models.utils.general.get_download"
    patch = mocker.patch(function, new=return_file)
    return patch


@pytest.fixture(autouse=True)
@pytest.mark.usefixtures("mocker")
def mock_upload(mocker, prepare_env):
    """Sending request to mock-minio."""
    function = "osp.models.utils.general.get_upload"
    patch = mocker.patch(function, new=store_file)
    return patch
