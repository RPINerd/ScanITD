""""""

from scanitd import __PACKAGE_NAME__, __version__


def test_package_name_is_defined() -> None:
    assert __PACKAGE_NAME__ == "ScanITD"


def test_version_is_semver_like() -> None:
    parts = __version__.split(".")
    assert len(parts) == 3
    assert all(part.isdigit() for part in parts)
