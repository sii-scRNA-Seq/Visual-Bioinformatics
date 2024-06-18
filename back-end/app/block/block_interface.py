from anndata import AnnData
from block.exception.missing_param_exception import MissingParametersException
import json


def adata_text(adata: AnnData) -> str:
    return f"Object with: {adata.n_obs:,} cells and {adata.n_vars:,} genes"


class Block:
    required_parameters = []

    def __init__(self):
        pass

    def validate_parameters(self, parameters: dict) -> None:
        missing_parameters = []
        for param in self.required_parameters:
            if param not in parameters:
                missing_parameters.append(param)
        if missing_parameters:
            raise MissingParametersException("Missing parameters: " + json.dumps(missing_parameters))

    def run(self, adata: AnnData, parameters: dict) -> (AnnData, dict):
        raise NotImplementedError()
