from anndata import AnnData
from block.block_interface import Block, adata_text
import scanpy as sc


class BasicFiltering(Block):

    """
    `BasicFiltering` subclass, which inherits from the `Block` superclass.
    """

    required_parameters = ["min_genes", "min_cells"]
    """The parameters required by a `BasicFiltering` block."""

    def __init__(self):
        """Initialise a `BasicFiltering` object."""
        super().__init__()

    def validate_parameters(self, parameters: dict) -> None:
        """
        Validate that all of the parameters in the `required_parameters` attribute are present in the parameters dictionary. The implementation is inherited from the `Block` superclass.

        Parameters:

            - `parameters`: A dictionary mapping block parameters to their values.
        """
        super(BasicFiltering, self).validate_parameters(parameters)

    def run(self, adata: AnnData, parameters: dict) -> (AnnData, dict):
        """
        Execute the code for a `BasicFiltering` block.

        Extract the values for `min_genes` and `min_cells` from the `parameters` dictionary, then filter cells and genes using Scanpy functions.

        Parameters:

            - `adata`: The AnnData for which the code should be executed.
            - `parameters`: A dictionary mapping parameter names to their values, which should be used while executing the code.

        Return:

            - The resulting AnnData after performing the block's behaviour.
            - A dictionary containing the results that will be seen by the user.
        """
        self.validate_parameters(parameters)
        min_genes = int(parameters["min_genes"])
        min_cells = int(parameters["min_cells"])

        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)

        message = {
            "text": adata_text(adata)
        }
        return adata, message
