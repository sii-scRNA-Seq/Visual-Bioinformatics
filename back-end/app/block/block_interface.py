from anndata import AnnData
from block.exception.missing_param_exception import MissingParametersException
import json


def adata_text(adata: AnnData) -> str:
    return f"Object with: {adata.n_obs:,} cells and {adata.n_vars:,} genes"


class Block:

    """
    Block superclass. Should only be instantiated through subclasses and never directly.

    Contains two methods:

        - validate_parameters() validates that all of the required parameters are present in the given dictionary and is provided concretely.
        - run() executes the code expected for a particular block and is provided abstractly.
    """

    required_parameters = []
    """The parameters required by a particular block."""


    def __init__(self):
        """Initialise a `Block` object."""
        pass

    def validate_parameters(self, parameters: dict) -> None:
        """
        Validate that all of the parameters in the required_parameters attribute are also present in the parameters dictionary. Raise an exception if this is not the case.

        Parameters:

            - `parameters`: a dictionary mapping block parameters to their values.

        Exceptions: `MissingParametersException`.
        """
        missing_parameters = []
        for param in self.required_parameters:
            if param not in parameters:
                missing_parameters.append(param)
        if missing_parameters:
            raise MissingParametersException("Missing parameters: " + json.dumps(missing_parameters))

    def run(self, adata: AnnData, parameters: dict) -> (AnnData, dict):
        """
        Execute the code for a particular block.

        Parameters:

            - `adata`: the AnnData for which the code should be executed.
            - `parameters`: a dictionary mapping parameter names to their values, which should be used while executing the code.

        Return:

            - The resulting AnnData after performing the block's behaviour.
            - A dictionary containing the results that will be seen by the user.
        """
        raise NotImplementedError()
