from anndata import AnnData
from block.block_interface import Block, adata_text
import scanpy as sc


class QCFiltering(Block):

    required_parameters = ["sample", "min_n_genes_by_counts", "max_n_genes_by_counts", "pct_counts_mt"]

    def __init__(self):
        super().__init__()

    def validate_parameters(self, parameters: dict) -> None:
        super(QCFiltering, self).validate_parameters(parameters)

    def run(self, adata: AnnData, dataset: str, parameters: dict) -> (AnnData, dict):
        self.validate_parameters(parameters)
        sample = str(parameters["sample"])
        min_n_genes_by_counts = float(parameters["min_n_genes_by_counts"])
        max_n_genes_by_counts = float(parameters["max_n_genes_by_counts"])
        pct_counts_mt = float(parameters["pct_counts_mt"])

        if dataset == "pbmc3k":
            adata.var["mt"] = adata.var_names.str.startswith("MT-")
            adata.obs["sample"] = "1"
        elif dataset == "pf_dogga":
            adata.var["mt"] = adata.var_names.str.contains("MIT")
            adata.obs["sample"] = adata.obs["day"]
        else:
            raise Exception("Selected dataset does not exist.")

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
