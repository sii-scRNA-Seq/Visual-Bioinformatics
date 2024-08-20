from anndata import AnnData
from block.block_interface import Block, adata_text
import scanpy as sc


class QCFiltering(Block):

    """
    `QCFiltering` subclass, which inherits from the `Block` superclass.
    """

    required_parameters = ["sample", "min_n_genes_by_counts", "max_n_genes_by_counts", "pct_counts_mt"]
    """The parameters required by a `QCFiltering` block."""

    def __init__(self):
        """Initialise a `QCFiltering` object."""
        super().__init__()

    def validate_parameters(self, parameters: dict) -> None:
        """
        Validates that all of the parameters in the `required_parameters` attribute are present in the parameters dictionary. The implementation is inherited from the `Block` superclass.

        Parameters:

            - `parameters`: A dictionary mapping block parameters to their values.
        """
        super(QCFiltering, self).validate_parameters(parameters)

    def run(self, adata: AnnData, dataset: str, parameters: dict) -> (AnnData, dict):
        """
        Execute the code for a `QCFiltering` block.

        Extracts the values for `sample`, `min_n_genes_by_counts`, `max_n_genes_by_counts` and `pct_counts_mt` from the `parameters` dictionary, calculates the QC metrics, then filters cells from the chosen sample based on these metrics.

        Parameters:

            - `adata`: The AnnData for which the code should be executed.
            - `dataset`: The user's selected dataset.
            - `parameters`: A dictionary mapping parameter names to their values, which should be used while executing the code.

        Return:

            - The resulting AnnData after performing the block's behaviour.
            - A dictionary containing the results that will be seen by the user.
        """
        self.validate_parameters(parameters)
        sample = str(parameters["sample"])
        min_n_genes_by_counts = float(parameters["min_n_genes_by_counts"])
        max_n_genes_by_counts = float(parameters["max_n_genes_by_counts"])
        pct_counts_mt = float(parameters["pct_counts_mt"])

        if dataset == "pbmc3k":
            adata.var["mt"] = adata.var_names.str.startswith("MT-")
            adata.obs["sample"] = adata.obs["sample"]
        elif dataset == "pf_dogga":
            adata.var["mt"] = adata.var_names.str.contains("MIT")
            adata.obs["sample"] = adata.obs["day"]
        else:
            raise Exception("Selected dataset does not exist.")

        if sample not in adata.obs["sample"].unique():
            raise Exception("Selected sample does not exist.")

        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        adata.obs["total_UMIs"] = adata.obs["total_counts"]
        adata.obs = adata.obs.drop("total_counts", axis=1)

        adata = adata[(adata.obs["sample"] != sample) | ((adata.obs["sample"] == sample) & (adata.obs["n_genes_by_counts"] > min_n_genes_by_counts))]
        adata = adata[(adata.obs["sample"] != sample) | ((adata.obs["sample"] == sample) & (adata.obs["n_genes_by_counts"] < max_n_genes_by_counts))]
        adata = adata[(adata.obs["sample"] != sample) | ((adata.obs["sample"] == sample) & (adata.obs["pct_counts_mt"] < pct_counts_mt))]

        message = {
            "text": adata_text(adata)
        }
        return adata, message
