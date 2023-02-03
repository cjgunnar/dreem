
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

import os
from os import path
import glob

mapping = {'.': 0, '(': 1, ')': 2}


def get_name(file_path):
    file_name = path.split(file_path)[1]
    if "reference" in file_name:
        name = path.splitext(file_name)[0]
        return name

    # example file name
    # 4_WT_plus_plus_S4_R1_001-ecoli_ryhB-05-85-K2_Cluster2_expUp_200_expDown_200.dot
    condition_parts, _, _, _, drawing_parts = file_name.split(sep='-')
    num_clusters, cluster_num, *_ = drawing_parts.split(sep='_')

    id_str, strain, DIP, DMS, *_ = condition_parts.split(sep="_")

    id_int = int(id_str.rstrip("bB"))

    if 1 <= id_int < 13:
        replicate = 1
    elif 13 <= id_int < 24:
        replicate = 2
    elif 24 <= id_int < 37:
        replicate = 3
    else:
        raise Exception(f"{id_int} is not mapped to a replicate")

    # return ' '.join([strain, "BR " + str(replicate), DIP, DMS, num_clusters, cluster_num])
    return ' '.join([strain, num_clusters, cluster_num])


def read_dots(files):
    header_lines = 2
    arrays = []
    names = []
    num_refs = 0
    for file in files:
        with open(file, encoding='utf-8') as f:
            for i, line in enumerate(f.readlines()):
                if i < header_lines:
                    continue
                if len(line) < 0 or line[0] not in '.()':
                    continue

                # 1-based
                structure_number = (i - header_lines) // 2 + 1

                #if structure_number > 1 and "reference" not in file:
                #    continue

                if "reference" in file:
                    num_refs += 1

                names.append(get_name(file) + " struct " + str(structure_number))
                # names.append(get_name(file))

                # read the (). string and convert to array of integers
                arrays.append(np.fromiter((mapping[c] for c in line if c in mapping), dtype='int8'))

    result = pd.DataFrame.from_dict({key: value for key, value in zip(names, arrays)})

    # rows=samples, columns=nucleotide number
    result = result.transpose()
    # sort sample names but the keep the reference (which will stay at the bottom)
    result[:-num_refs] = result[:-num_refs].sort_index()

    return result


def main(files, reference_name, output_fig_name):
    # import data
    df = read_dots(files)

    # create figure
    fig, ax = plt.subplots(figsize=(10, 6.4))

    cmap = sns.color_palette(["white", "blue", "green"])
    ax = sns.heatmap(df, cbar=False, cmap=cmap, xticklabels=5, yticklabels=True, ax=ax)

    for i in range(df.shape[0] + 1):
        ax.axhline(i, color='white', lw=2)

    # create a legend
    # I've opted to make a diagram explaining it instead for the poster
    # single = patches.Patch(color="white", label="Single")
    # double_up = patches.Patch(color="blue", label="Double Rising")
    # double_down = patches.Patch(color="green", label="Double Falling")
    # plt.legend(handles=[single, double_up, double_down])

    ax.set_title(reference_name)
    ax.set_ylabel("strain DIP DMS clusters cluster# structure#")
    ax.set_xlabel("Nucleotide")

    # output figure
    fig.tight_layout()
    fig.savefig(output_fig_name, dpi=300)


def snakemake_run(snakemake):
    directories = snakemake.input
    reference_dot = snakemake.input.reference_dot
    output_file = snakemake.output[0]
    reference_name = snakemake.wildcards.ref

    files = []
    for directory in directories:
        dot_files = glob.glob(f"{directory}/*-*/K_*/*-*-*-*-K*_Cluster*_expUp_200_expDown_200.dot")
        files += dot_files

    files += [reference_dot]

    output_dir = path.split(output_file)[0]
    if not path.exists(output_dir):
        os.makedirs(output_dir)

    main(files, reference_name, output_file)


if __name__ == "__main__":
    snakemake_run(snakemake)
