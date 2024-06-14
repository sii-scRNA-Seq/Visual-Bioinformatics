import pytest
from basic_filtering import BasicFiltering
from exception.missing_param_exception import MissingParametersException


def test_param_validation_requires_min_genes():
    block = BasicFiltering()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "min_cells": "25"
        })
    assert e_info.value.message == 'Missing parameters: ["min_genes"]'


def test_param_validation_requires_min_cells():
    block = BasicFiltering()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "min_genes": "25"
        })
    assert e_info.value.message == 'Missing parameters: ["min_cells"]'

# TODO add tests for running block
# TODO remove now superfluous tests from e2e
