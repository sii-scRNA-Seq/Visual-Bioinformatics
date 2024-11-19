from anndata import AnnData
from block.block_interface import Block
import codecs
import io
from matplotlib import pyplot as plt
import scanpy as sc

from block.exception.invalid_param_exception import InvalidParametersException


class PlotDimensionReduction(Block):
    """
    `PlotReducedDimension` subclass, which inherits from the `Block` superclass.
    """

    required_parameters = ["reduction", "group_by"]
    """The parameters required by a `PlotDimensionReduction` block."""

    def __init__(self):
        """Initialise a `PlotReducedDimension` object."""
        super().__init__()

    def validate_parameters(self, parameters: dict) -> None:
        """
        Validate that all of the parameters in the `required_parameters` attribute are present in the parameters dictionary. The implementation is inherited from the `Block` superclass.

        Parameters:

            - `parameters`: A dictionary mapping block parameters to their values.
        """
        super(PlotDimensionReduction, self).validate_parameters(parameters)

    def run(self, adata: AnnData, dataset: str, parameters: dict) -> (AnnData, dict):
        """
        Execute the code for a `PlotDimensionReduction` block.

        Plot a precomputed UMAP or, integrated, or PCA embedding.

        Parameters:

            - `adata`: The AnnData for which the code should be executed.
            - `dataset`: The user's selected dataset.
            - `parameters`: A dictionary mapping parameter names to their values, which should be used while executing the code.

        Return:

            - The resulting AnnData after performing the block's behaviour.
            - A dictionary containing the results that will be seen by the user.
        """
        self.validate_parameters(parameters)

        reduction = str(parameters["reduction"])
        color_by = str(parameters["group_by"])

        color = None
        if color_by in list(adata.obs):
            color = color_by

        if reduction == "PCA":
            if "X_pca" in adata.obsm.keys():
                sc.pl.pca(adata, color=color, show=False)
            else:
                raise InvalidParametersException("PCA not run yet")
        elif reduction == "UMAP":
            if "X_umap" in adata.obsm.keys():
                sc.pl.umap(adata, color=color, show=False)
            else:
                raise InvalidParametersException("UMAP not run yet")
        elif reduction == "TSNE":
            if "X_tsne" not in adata.obsm.keys():
                sc.tl.tsne(adata)
            sc.pl.tsne(adata, color=color, show=False)
        else:
            raise InvalidParametersException(f"{reduction} not supported")

        image_stream = io.BytesIO()
        plt.savefig(image_stream, format="png")
        image_stream.seek(0)
        image = {
            "image": str(codecs.encode(image_stream.read(), "base64")),
            "alt_text": f"A scatter plot of the data, reduced by {reduction}."
        }
        message = {
            "image": image
        }
        return adata, message
