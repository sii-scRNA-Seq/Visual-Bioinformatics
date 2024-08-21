from anndata import AnnData
from block.exception.missing_param_exception import MissingParametersException
import json


def adata_text(adata: AnnData) -> str:
    """
    Generate text to describe an AnnData.

    Parameters:

        - `adata`: The AnnData from which text should be generated.

    Return:

        - Text describing the given AnnData.
    """
    return f"Object with: {adata.n_obs:,} cells and {adata.n_vars:,} genes"


class Block:

    """
    `Block` superclass. Should only be instantiated through subclasses and never directly.

    A block represents a stage in the single cell pipeline, responsible for receiving and validating parameters, then executing that stage of the pipeline.
    """

    required_parameters = []
    """The parameters required by a particular block, which will be used to validate that these parameters exist."""

    def __init__(self):
        """Initialise a `Block` object."""
        pass

    def validate_parameters(self, parameters: dict) -> None:
        """
        Validate that all of the parameters in the `required_parameters` attribute are present in the `parameters` dictionary. Raise an exception if this is not the case.

        Parameters:

            - `parameters`: A dictionary mapping block parameters to their values.

        Exceptions: `MissingParametersException`.
        """
        missing_parameters = []
        for param in self.required_parameters:
            if param not in parameters:
                missing_parameters.append(param)
        if missing_parameters:
            raise MissingParametersException("Missing parameters: " + json.dumps(missing_parameters))

    def run(self, adata: AnnData, dataset: str, parameters: dict) -> (AnnData, dict):
        """
        Execute the code for a particular block.

        Parameters:

            - `adata`: The AnnData for which the code should be executed.
            - `dataset`: The user's selected dataset.
            - `parameters`: A dictionary mapping parameter names to their values, which should be used while executing the code.

        Return:

            - The resulting AnnData after performing the block's behaviour.
            - A dictionary containing the results that will be seen by the user.

        Exceptions: `NotImplementedError`.
        """
        raise NotImplementedError()
