"""Helper functions for OSP-utilities."""
import io
import json
import os
import tempfile
import warnings
from typing import Any, Dict, List, TextIO, Union
from uuid import uuid4

import requests
from minio import Minio
from minio.error import MinioException
from urllib3.response import HTTPResponse

import osp.dictionaries as cases
from osp.models.settings.general import ModelSettings

CASES = os.path.dirname(cases.__file__)


def get_minio() -> Minio:
    """Instantiate MinIO-client without `Depends`."""
    settings = ModelSettings()
    return Minio(
        settings.minio_endpoint,
        access_key=settings.minio_user.get_secret_value(),
        secret_key=settings.minio_password.get_secret_value(),
        secure=False,
    )


def _get_upload(filepath: str, uuid: str, minio_client: Minio) -> str:
    """Helper function for `depends_upload` and `get_upload`."""
    suffix = os.path.splitext(filepath)[-1]
    # Generate a unique UUID as the object key
    object_key = uuid or str(uuid4())

    # Upload the file to MinIO
    try:
        if not minio_client.bucket_exists(object_key):
            minio_client.make_bucket(object_key)
        minio_client.fput_object(
            object_key,
            object_key,
            filepath,
            metadata={"suffix": suffix},
        )
    except Exception as err:
        raise MinioException(err.args) from err
    return object_key


def _get_download(uuid: str, minio_client: Minio) -> HTTPResponse:
    """Helper function for `depends_download` and `get_download`."""
    file_data = io.BytesIO()
    if not minio_client.bucket_exists(uuid):
        raise MinioException("Bucket does not exist.")
    try:
        response = minio_client.get_object(
            uuid,
            uuid,
            file_data,
        )
    except Exception as err:
        raise MinioException(err) from err
    return response


def get_upload(file: Union[str, TextIO], uuid: str = None) -> str:
    """Upload file with minio client and return upload-id without `Depends`."""
    if hasattr(file, "filename"):
        file = file.filename
    elif hasattr(file, "name"):
        file = file.name
    minio_client = get_minio()
    return _get_upload(file, uuid, minio_client)


def get_download(uuid: str, as_file=False) -> Union[HTTPResponse, str]:
    """Return `io.BytesIO` from uuid through minio client without `Depends`."""
    minio_client = get_minio()
    response = _get_download(uuid, minio_client)
    if as_file:
        suffix = response.headers.get("x-amz-meta-suffix")
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as temp:
            temp.write(response.data)
            response = temp.name
    return response


def _get_example_json(example_json: str, standard_data: List) -> Dict[str, Any]:
    try:
        for uuid, filename in standard_data:
            try:
                get_download(uuid)
            except Exception:
                filepath = os.path.join(CASES, "example_data", filename)
                with open(filepath, "rb") as file:
                    get_upload(file, uuid=uuid)
    except Exception as error:
        warnings.warn(f"Error while uploading the example file to MinIO: {error.args}")
    return _get_example("ams", example_json)


def _get_example(folder: str, file: str) -> Dict[str, Any]:
    filepath = os.path.join(CASES, folder, file)
    with open(filepath, "r+") as handler:
        return json.load(handler)


def _download_file(url) -> str:
    settings = ModelSettings()
    response = requests.get(url, timeout=settings.download_timeout)
    return response.text
