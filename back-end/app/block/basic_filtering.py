from anndata import AnnData
from block.block_interface import Block, adata_text
import scanpy as sc


class BasicFiltering(Block):

    required_parameters = ["min_genes", "min_cells"]

    def __init__(self):
        super().__init__()

    def validate_parameters(self, parameters: dict) -> None:
        super(BasicFiltering, self).validate_parameters(parameters)

    def run(self, adata: AnnData, parameters: dict) -> (AnnData, dict):
        self.validate_parameters(parameters)
        min_genes = int(parameters["min_genes"])
        min_cells = int(parameters["min_cells"])

        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)

        message = {
            "text": adata_text(adata)
        }
        return adata, message
