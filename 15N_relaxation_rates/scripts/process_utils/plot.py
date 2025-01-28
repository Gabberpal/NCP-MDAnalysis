import os
import matplotlib.pyplot as plt


def get_autocorr_graph_label(fit_line):
    """
    Generates a formatted label for an autocorrelation graph based on fitted parameters.

    This function takes a pandas Series (`fit_line`) containing fitted parameters (amplitudes and relaxation times)
    and creates a human-readable label for a graph. The label includes the amplitude (A) and relaxation time (τ) values
    in a formatted string.

    Args:
        fit_line (pd.Series): A pandas Series containing fitted parameters. The Series should include keys
                              for amplitudes (ending with '-a') and relaxation times (ending with '-tau').

    Returns:
        str: A formatted string representing the graph label. Each line contains an amplitude and its corresponding
             relaxation time in the format:
             "A = <value> ; τ = <value>"
    """
    # Extract amplitude and tau values from the Series
    amplitude = fit_line.filter(like='-a') # Filter keys containing '-a' (amplitudes)
    tau = fit_line.filter(like='-tau') # Filter keys containing '-tau' (relaxation times)

    # Create a list of formatted strings for each amplitude-tau pair
    union_a_tau = ["{a_label:2s} = {a_value:5.3e} ; {tau_label:3s} = {tau_value: 8.3e}".format(
        a_label=a_label,
        a_value=fit_line[a_label],
        tau_label=tau_label,
        tau_value=fit_line[tau_label])
        for a_label, tau_label in zip(amplitude.index.tolist(), tau.index.tolist())
    ]

    # Join the formatted strings into a single multi-line label
    graph_label = "\n".join(union_a_tau)
    return graph_label


def add_relpath_to_top_corner(figure: plt.Figure):
    """
    Adds a relative path of the current working directory (relative to a specified base directory)
    to the top-right corner of the given matplotlib figure.

    Args:
        figure (plt.Figure): The matplotlib figure to which the text will be added.

    The function creates a transparent axis covering the entire figure and adds a text annotation
    in the top-right corner with the relative path of the current directory.
    """
    # Create a transparent axis that covers the entire figure
    big_axis = figure.add_axes([0, 0, 1, 1], facecolor=(1, 1, 1, 0))

    # Calculate the relative path of the current directory relative to the specified base directory
    relpath = os.path.relpath(os.path.abspath(os.curdir), os.path.expanduser("~/bioinf/handling"))

    # Add the relative path as text to the top-right corner of the figure
    big_axis.text(
        0.99,  # x-coordinate (right-aligned)
        0.99,  # y-coordinate (top-aligned)
        relpath,  # Text to display (relative path)
        color="#CCCCCC",  # Light gray color for the text
        horizontalalignment='right',  # Align text to the right
        verticalalignment='top'  # Align text to the top
    )


def settings_plot(graph_label):
    left, width = .40, .54
    bottom, height = .40, .54
    right = left + width
    top = bottom + height
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.text(right, top, graph_label,
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes,
            multialignment='left',
            bbox={'facecolor': 'moccasin', 'alpha': 0.5, 'pad': 6})
    ax.set_xlabel('time, ns', fontsize=13)
    ax.set_ylabel('C(t)', fontsize=13)
    return fig, ax


def set_axis_parameters(xname, yname, title_tag=""):
    fig = plt.figure(figsize=(20, 12))
    ax = fig.add_subplot(111)
    ax.set_xlabel('{xname}'.format(xname=xname), fontsize=22)
    ax.set_ylabel('{yname}'.format(yname=yname), fontsize=22)
    ax.set_title('{title_tag}'.format(title_tag=title_tag), fontsize=25, loc='center')
    plt.xticks(fontsize=19)
    plt.yticks(fontsize=19)
    ax.grid(True)
    return fig, ax
