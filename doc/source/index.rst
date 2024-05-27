:html_theme.sidebar_secondary.remove: true

.. template taken from Pandas

.. module:: scipy

*******************
SciPy documentation
*******************

**Date**: |today| **Version**: |version|

**Download documentation**: https://docs.scipy.org/doc/

**Useful links**:
`Install <https://scipy.org/install/>`__ |
`Source Repository <https://github.com/scipy/scipy>`__ |
`Issues & Ideas <https://github.com/scipy/scipy/issues>`__ |
`Q&A Support <https://stackoverflow.com/questions/tagged/scipy>`__ |
`Mailing List <https://mail.python.org/mailman3/lists/scipy-dev.python.org/>`__

**SciPy** (pronounced "Sigh Pie") is an open-source software for mathematics,
science, and engineering.

.. grid:: 1 1 2 2
    :gutter: 2 3 4 4

    .. grid-item-card::
        :img-top: _static/index_user_guide.svg
        :text-align: center

        **User guide**
        ^^^

        The user guide provides in-depth information on the
        key concepts of SciPy with useful background information and explanation.

        +++

        .. button-ref:: user_guide
            :color: secondary
            :click-parent:

            To the user guide

    .. grid-item-card::
        :img-top: _static/index_api.svg
        :text-align: center

        **API reference**
        ^^^

        The reference guide contains a detailed description of
        the SciPy API. The reference describes how the methods work and which parameters can
        be used. It assumes that you have an understanding of the key concepts.

        +++

        .. button-ref:: scipy-api
            :color: secondary
            :click-parent:

            To the reference guide

    .. grid-item-card::
        :img-top: _static/index_getting_started.svg
        :text-align: center

        **Building from source**
        ^^^

        Want to build from source rather than use a Python distribution or
        pre-built SciPy binary? This guide will describe how to set up your
        build environment, and how to build *SciPy* itself, including the many
        options for customizing that build.

        +++

        .. button-ref:: building-from-source
            :color: secondary
            :click-parent:

            To the build guide

    .. grid-item-card::
        :img-top: _static/index_contribute.svg
        :text-align: center

        **Developer guide**
        ^^^

        Saw a typo in the documentation? Want to improve
        existing functionalities? The contributing guidelines will guide
        you through the process of improving SciPy.

        +++

        .. button-ref:: scipy-development
            :color: secondary
            :click-parent:

            To the development guide

.. toctree::
   :maxdepth: 1
   :hidden:

   Installing <https://scipy.org/install/>
   User Guide <tutorial/index>
   API reference <reference/index>
   Building from source <building/index>
   Development <dev/index>
   Release notes <release>
