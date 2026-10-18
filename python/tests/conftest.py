"""Optional tests requiring non-public input data are opt-in."""

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--private-data", action="store_true", help="Run optional private-data tests"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--private-data"):
        return
    skip = pytest.mark.skip(reason="private data: opt in with --private-data")
    for item in items:
        if "private_data" in item.keywords:
            item.add_marker(skip)
