from IPython.core.pylabtools import print_figure
import matplotlib.pyplot as plt
import numpy as np


def _repr_png_(distribution):

    name = distribution.dist.name
    title = (f"{name}(args={distribution.args}, kwargs={distribution.kwds}),"
             " N=1000")

    fig, ax = plt.subplots()

    # PDF
    x = np.linspace(distribution.ppf(0.01),
                    distribution.ppf(0.99), 100)
    pdf = distribution.pdf(x)
    ax.plot(x, pdf, '-', lw=5, alpha=0.6)

    # Empirical PDF
    sample = distribution.rvs(size=1000)
    ax.hist(sample, density=True, histtype="stepfilled", alpha=0.2)
    delta = np.max(pdf) * 5e-2
    ax.plot(sample[:100], -delta - delta * np.random.random(100), "+k")

    ax.set_title(title)
    ax.set_ylabel(r"$f$")
    ax.set_xlabel(r"$x$")

    data = print_figure(fig, "png")

    # We MUST close the figure, otherwise IPython's display machinery
    # will pick it up and send it as output, resulting in a double display
    plt.close(fig)
    return data


def load_ipython_extension(ipython):
    png_f = ipython.display_formatter.formatters["image/png"]
    png_f.for_type_by_name("scipy.stats._distn_infrastructure",
                           "rv_frozen", _repr_png_)
