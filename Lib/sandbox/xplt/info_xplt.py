"""
Plotting routines based on Gist 
===============================

  All Gist functions are also available as xplt.<command>
  xplt.ghelp('plsys')

    histogram -- Plot a histogram.
    barplot   -- Construct a barplot.
    errorbars -- Draw connected points with (y-only) errorbars.
    legend    -- Construct and place a legend.
    arrow     -- Draw an arrow.
    plot      -- Plot curves.
    logxy     -- Change scale to logarithmic
    surf      -- Plot surfaces (no axes currently plotted).
    xlabel    -- Place a label on the x-axis.
    ylabel    -- Place a label on the y-axis.
    title     -- Place a title above the plot.
    hold      -- Draw subsequent plots over the current plot.
    matplot   -- Plot many curves in a matrix against a single x-axis.
    addbox    -- Add a box to the current plot.
    imagesc   -- Draw an image.
    imagesc_cb -- Draw an image with a colorbar.
    movie     -- Play a sequence of images as a movie.
    figure    -- Create a new figure.
    full_page -- Create a full_page window.
    subplot   -- Draw a sub-divided plot.
    plotframe -- Change the plot system on a multi-plot page.
    twoplane  -- Create a plot showing two orthogonal planes of a
                 three-dimensional array.

Xplt is basically a wrapper around the lower-level gist plotting library
used with Yorick.  Simpler, higher-level commands are provided in this
module, but all of the low-level commands of gist are still available under
the xplt namespace.

See

  http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/python/pygist_html/pygist.html

for more information on the gist package and it's commands.
"""

postpone_import = 1
