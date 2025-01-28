import os
import matplotlib.pyplot as plt


def get_autocorr_graph_label(fit_line):
    amplitude = fit_line.filter(like='-a')
    tau = fit_line.filter(like='-tau')
    union_a_tau = ["{a_label:2s} = {a_value:5.3e} ; {tau_label:3s} = {tau_value: 8.3e}".format(
        a_label=a_label,
        a_value=fit_line[a_label],
        tau_label=tau_label,
        tau_value=fit_line[tau_label])
        for a_label, tau_label in zip(amplitude.index.tolist(), tau.index.tolist())
    ]
    graph_label = "\n".join(union_a_tau)
    return graph_label


def add_relpath_to_top_corner(figure: plt.Figure):
    big_axis = figure.add_axes([0, 0, 1, 1], facecolor=(1, 1, 1, 0))

    big_axis.text(0.99, 0.99,
                  os.path.relpath(os.path.abspath(os.curdir), os.path.expanduser("~/bioinf/handling")),
                  color="#CCCCCC",
                  horizontalalignment='right',
                  verticalalignment='top')


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
