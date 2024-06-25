from anndata import AnnData
from block.block_interface import Block, adata_text
from dataset_info import dataset_info
import scanpy as sc


class Integration(Block):

    required_parameters = ["observation"]

    def __init__(self):
        super().__init__()

    def validate_parameters(self, parameters: dict) -> None:
        super(Integration, self).validate_parameters(parameters)

    def run(self, adata: AnnData, dataset: str, parameters: dict) -> (AnnData, dict):
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
