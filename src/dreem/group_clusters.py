
import os
import glob
from shutil import copytree


def group_clusters(range_dirs, output):
    """from a list of directories with cluster results, move the clusters together"""

    for range_dir in range_dirs:
        best_run_dirs = glob.glob(os.path.join(range_dir, "K_*/run_*-best"), recursive=True)
        for best_run_dir in best_run_dirs:
            # brittle but effective way at getting K_2 in the
            # path to the best_run_directory
            *_, cluster_count, _ = best_run_dir.split(sep="/")
            *_, start, end = os.path.basename(range_dir).split("-")
            destination = os.path.join(output, f"{start}-{end}", cluster_count)
            # print(f"copy {best_run_dir} to {destination}")
            copytree(best_run_dir, destination, dirs_exist_ok=True)


# noinspection PyUnresolvedReferences
group_clusters(snakemake.input, snakemake.output[0])
