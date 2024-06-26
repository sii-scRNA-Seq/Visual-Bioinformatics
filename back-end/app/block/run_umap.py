from anndata import AnnData
from block.block_interface import Block
import codecs
import io
from joblib import parallel_backend
from matplotlib import pyplot as plt
import scanpy as sc
from threadpoolctl import threadpool_limits


class RunUMAP(Block):

    required_parameters = ["n_neighbors", "n_pcs"]

    def __init__(self):
        super().__init__()

    def validate_parameters(self, parameters: dict) -> None:
        super(RunUMAP, self).validate_parameters(parameters)

    def run(self, adata: AnnData, parameters: dict) -> (AnnData, dict):
        self.validate_parameters(parameters)

        n_neighbors = int(parameters["n_neighbors"])
        n_pcs = int(parameters["n_pcs"])

        with parallel_backend("threading", n_jobs=1):
            with threadpool_limits(limits=1, user_api="blas"):
                sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
                sc.tl.umap(adata)
                sc.tl.leiden(adata)
                sc.pl.umap(adata, color=["leiden"], legend_loc="on data")

        image_stream = io.BytesIO()
        plt.savefig(image_stream, format="png")
        image_stream.seek(0)
        message = {
            "img": str(codecs.encode(image_stream.read(), "base64")),
            "alttext": "A UMAP of Leiden clusters using the principal components generated by a Run UMAP block"
        }
        return adata, message
