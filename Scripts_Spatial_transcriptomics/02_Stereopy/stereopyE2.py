#!/usr/bin/env python3

import argparse
import warnings
import stereo as st
import os
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


def save_plot(plot_dir, filename, dpi=300):
    """Save current matplotlib figure."""
    plt.savefig(
        os.path.join(plot_dir, filename),
        dpi=dpi,
        bbox_inches="tight"
    )
    plt.close()


def main():

    parser = argparse.ArgumentParser(
        description="StereoPy pipeline (Snakemake single-bin execution)"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input tissue.gef file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output h5ad file"
    )

    parser.add_argument(
        "--bin_size",
        type=int,
        required=True,
        help="Bin size"
    )

    parser.add_argument(
        "--npcs",
        type=int,
        default=30,
        help="Number of principal components"
    )

    args = parser.parse_args()

    print(f"Stereo version : {st.__version__}")
    print(f"Input          : {args.input}")
    print(f"Bin size       : {args.bin_size}")

    ####################################################
    # Output folders
    ####################################################

    outdir = os.path.dirname(args.output)

    plot_dir = os.path.join(outdir, "plots")
    marker_dir = os.path.join(outdir, "marker_genes")

    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(marker_dir, exist_ok=True)

    ####################################################
    # Read data
    ####################################################

    data = st.io.read_gef(
        file_path=args.input,
        bin_type="bins",
        bin_size=args.bin_size,
        is_sparse=True
    )

    print(data)

    ####################################################
    # QC
    ####################################################

    data.tl.cal_qc()

    # Example (depends on StereoPy version)
    # data.plt.violin()
    # save_plot(plot_dir, "QC.png")

    data.tl.filter_cells(
        min_gene=20,
        min_n_genes_by_counts=3,
        pct_counts_mt=15,
        inplace=True
    )

    data.tl.raw_checkpoint()

    ####################################################
    # Normalization
    ####################################################

    data.tl.normalize_total(target_sum=10000)
    data.tl.log1p()

    ####################################################
    # Highly variable genes
    ####################################################

    data.tl.highly_variable_genes(
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_top_genes=2000,
        res_key="highly_variable_genes"
    )

    ####################################################
    # Scaling
    ####################################################

    data.tl.scale(
        max_value=10,
        zero_center=True
    )

    ####################################################
    # PCA
    ####################################################

    data.tl.pca(
        use_highly_genes=True,
        n_pcs=args.npcs,
        res_key="pca"
    )

    ####################################################
    # Neighbors
    ####################################################

    data.tl.neighbors(
        pca_res_key="pca",
        n_pcs=args.npcs,
        res_key="neighbors"
    )

    data.tl.spatial_neighbors(
        neighbors_res_key="neighbors",
        res_key="spatial_neighbors"
    )

    ####################################################
    # UMAP
    ####################################################

    data.tl.umap(
        pca_res_key="pca",
        neighbors_res_key="neighbors",
        res_key="umap"
    )

    ####################################################
    # Clustering
    ####################################################

    data.tl.leiden(
        neighbors_res_key="neighbors",
        res_key="leiden"
    )

    data.tl.leiden(
        neighbors_res_key="spatial_neighbors",
        res_key="spatial_leiden"
    )

    ####################################################
    # Marker genes
    ####################################################

    data.tl.find_marker_genes(
        cluster_res_key="leiden",
        method="t_test",
        use_highly_genes=False,
        use_raw=True
    )

    markers = data.tl.result["marker_genes"]

    print("Marker gene object:", type(markers))

    if isinstance(markers, pd.DataFrame):

        markers.to_csv(
            os.path.join(marker_dir, "all_marker_genes.csv"),
            index=False
        )

        if "group" in markers.columns and "scores" in markers.columns:

            top10 = (
                markers
                .sort_values(["group", "scores"],
                             ascending=[True, False])
                .groupby("group")
                .head(10)
            )

            top10.to_csv(
                os.path.join(marker_dir,
                             "top10_marker_genes.csv"),
                index=False
            )

    elif isinstance(markers, dict):

        for key, value in markers.items():

            if isinstance(value, pd.DataFrame):

                value.to_csv(
                    os.path.join(marker_dir,
                                 f"{key}.csv"),
                    index=False
                )

    ####################################################
    # Save cluster assignments
    ####################################################

    try:

        cluster_df = pd.DataFrame(index=data.cells.cell_name)

        cluster_df["leiden"] = data.cells["leiden"]
        cluster_df["spatial_leiden"] = data.cells["spatial_leiden"]

        cluster_df.to_csv(
            os.path.join(outdir, "clusters.csv")
        )

    except Exception as e:

        print("Could not save cluster table:", e)

    ####################################################
    # Export h5ad
    ####################################################

    st.io.stereo_to_anndata(
        data,
        flavor="seurat",
        output=args.output
    )

    print("\nPipeline completed successfully.")
    print("Output:", args.output)


if __name__ == "__main__":
    main()
