from anndata import AnnData
from block.block_interface import Block, adata_text
from dataset_info import dataset_info
import scanpy as sc


class Integration(Block):

    """
    `Integration` subclass, which inherits from the `Block` superclass.
    """

    required_parameters = ["observation"]
    """The parameters required by an `Integration` block."""

    def __init__(self):
        """Initialise an `Integration` object."""
        super().__init__()

    def validate_parameters(self, parameters: dict) -> None:
        """
        Validates that all of the parameters in the `required_parameters` attribute are present in the parameters dictionary. The implementation is inherited from the `Block` superclass.

        Parameters:

            - `parameters`: A dictionary mapping block parameters to their values.
        """
        super(Integration, self).validate_parameters(parameters)

    def run(self, adata: AnnData, dataset: str, parameters: dict) -> (AnnData, dict):
        """
        Execute the code for an `Integration` block.

        Extracts the value for `observation` from the `parameters` dictionary, then integrates samples based on this observation.

        Parameters:

            - `adata`: The AnnData for which the code should be executed.
            - `dataset`: The user's selected dataset.
            - `parameters`: A dictionary mapping parameter names to their values, which should be used while executing the code.

        Return:

            - The resulting AnnData after performing the block's behaviour.
            - A dictionary containing the results that will be seen by the user.
        """
        self.validate_parameters(parameters)
        observation = str(parameters["observation"])

        if dataset not in [d["key"] for d in dataset_info]:
            raise Exception("Selected dataset does not exist.")
        for d in dataset_info:
            if d["key"] == dataset and observation not in d["integration_obs"]:
                raise Exception("Selected observation does not exist.")

        sc.external.pp.harmony_integrate(adata, observation)
        adata.obsm["X_pca"] = adata.obsm["X_pca_harmony"]

        message = {
            "text": adata_text(adata)
        }
        return adata, message
