import os
import matplotlib.pyplot as plt
import pandas as pd


def plot_and_save_rmsd(path_to_rmsd_dir: str,
                       output_directory: str = ".",
                       output_name: str = "rmsd.png"
                       ) -> None:
    # Указываем только один файл
    rmsd_column_names = {
        "rmsd_protein": "sec. struct. Cα",
        "rmsd_dna": "DNA",
        "rmsd_all": "sec. struct. Cα & DNA",
        "rmsd_dna_inner": "DNA inner turn",
        "rmsd_dna_outer": "DNA outer turn"
    }

    SMALL_SIZE = 16
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 20

    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)  # legend fontsize

    colors = ["blue", "red", "black", "green", "magenta"]
    fig, axs = plt.subplots(1, 2, figsize=(41 / 2.54, 12 / 2.54), dpi=300)

    axis_limits = None
    for ind, (key, label) in enumerate(rmsd_column_names.items()):
        df_rmsd = pd.read_csv(os.path.join(path_to_rmsd_dir, 'rmsd.csv'))  # Use the provided path
        axs[ind // 3].plot(df_rmsd["time_ns"], df_rmsd[key],
                           color=colors[ind], linewidth=1.0, label=label)
        if axis_limits is None:
            axis_limits = [0, df_rmsd["time_ns"].values[-1] + 10, 0, 14]

    for index in [0, 1]:
        axs[index].set(xlabel="time, ns", ylabel=r"RMSD, ${\rm\AA}$")
        axs[index].set_xlim(axis_limits[:2])
        axs[index].set_ylim(axis_limits[2:])

    ax_right_1 = axs[0].twinx()
    ax_right_1.set_yticks(axs[0].get_yticks())
    ax_right_1.set_ylim(axs[0].get_ylim())
    ax_right_1.set_yticklabels([])
    ax_right_2 = axs[1].twinx()
    ax_right_2.set_yticks(axs[1].get_yticks())
    ax_right_2.set_ylim(axs[1].get_ylim())
    ax_right_2.set_yticklabels([])
    leg_1 = axs[0].legend(loc="upper left", bbox_to_anchor=(0.305, 0.97))
    leg_2 = axs[1].legend(loc="upper left", bbox_to_anchor=(0.45, 0.97))

    for line in leg_1.get_lines():
        line.set_linewidth(4.0)

    for line in leg_2.get_lines():
        line.set_linewidth(4.0)

    os.makedirs(output_directory, exist_ok=True)
    plt.savefig(os.path.join(output_directory, output_name), bbox_inches="tight")
    plt.close()
