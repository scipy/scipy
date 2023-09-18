from typing import NamedTuple, Optional


class URL(NamedTuple):
    protocol: Optional[str]
    slashes: bool
    auth: Optional[str]
    port: Optional[str]
    hostname: Optional[str]
    hash: Optional[str]  # noqa: A003
    search: Optional[str]
    pathname: Optional[str]
