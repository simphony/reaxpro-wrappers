"""Model settings for OSP Pydantic Models"""
from pydantic import BaseSettings, Field, SecretStr


class ModelSettings(BaseSettings):
    """Main configuration for based on redis-settings"""

    minio_user: SecretStr = Field("rootname", description="User to MinIO instance.")
    minio_password: SecretStr = Field(
        "foobar123", description="Password to MinIO instance."
    )
    minio_endpoint: str = Field(
        "localhost:9000",
        description="Resolvable endpoint to contact MinIO instance.",
    )
    download_timeout: int = Field(
        10, description="Timeout when a file needs to be downloaded via http."
    )

    class Config:
        """Pydantic Config for ReaxProSettings."""

        env_prefix = "REAXPRO_"
