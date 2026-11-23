=================================
Contributing to the documentation
=================================

The documentation is built with `Sphinx <https://www.sphinx-doc.org/>`__, the `PyData Sphinx theme <https://pydata-sphinx-theme.readthedocs.io/>`__ and the `Sphinx-Gallery <https://sphinx-gallery.github.io/stable/index.html>`__ extension.

Layout of the sources
=====================

- ``doc/`` contains the Sphinx configuration (``conf.py``), the landing page (``index.rst``) and the other pages, written in reStructuredText. A new page must be added to the ``toctree`` at the end of ``index.rst`` to appear in the navigation bar.
- ``doc/_static/`` contains the logo, the favicon, the gallery thumbnails and ``switcher.json``, which lists the ``stable`` and ``dev`` versions for the version switcher.
- ``examples/`` in the root of the repository contains the sources of the *Tutorial and examples* gallery. ``examples/plot_data/`` holds the data files read by the executed examples.

The gallery pages are generated in ``doc/tutorial_and_examples/`` during the build. That folder is ignored by git and must not be edited.

Building locally
================

Install the requirements in a Python environment:

.. code:: bash

   python -m pip install sphinx pydata-sphinx-theme sphinx-gallery matplotlib numpy

``matplotlib`` and ``numpy`` are needed to execute the ``plot_*`` examples. Then build the HTML pages from the ``doc/`` folder:

.. code:: bash

   make html

and open ``_build/html/index.html`` in a browser. Check that the build prints no warnings.

Adding a tutorial or an example
===============================

Each example is a Python file in ``examples/``. ``examples/GALLERY_HEADER.rst`` is the landing page of the gallery, followed by the list of examples.
The file starts with a docstring in reStructuredText that holds the title and the text of the page; see the `Sphinx-Gallery syntax <https://sphinx-gallery.github.io/stable/syntax.html#python-script-syntax>`__ for how to alternate text and code blocks.
There are two kinds of files:

- files whose name starts with ``plot_``, like ``plot_direct_detection.py``, are executed during the build, and their output and figures are included in the page;
- all other files, like ``getting_started.py``, are rendered without being executed.

Every file in ``examples/`` becomes a gallery page, so keep standalone tools elsewhere (for example in ``tests/tools/``).
To choose the gallery thumbnail, add a PNG to ``doc/_static/thumbnail/`` and reference it at the end of the example:

.. code:: python

   # sphinx_gallery_thumbnail_path = '_static/thumbnail/my_example_thumb.png'

Publishing
==========

The workflow ``.github/workflows/build-and-deploy.yml`` builds the documentation on every push to ``dev`` or ``master`` that changes ``doc/`` or ``examples/``, and publishes it on the ``gh-pages`` branch: ``master`` goes to the `stable <https://maddmhep.github.io/maddm/stable/>`__ version and ``dev`` to the `dev <https://maddmhep.github.io/maddm/dev/>`__ version.
No manual copy to ``gh-pages`` is needed.
