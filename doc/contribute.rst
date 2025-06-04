=================================
Contributing to the documentation
=================================

The documentation has been created with `Sphinx <https://www.sphinx-doc.org/>`__.

How this documentation was created
==================================

Install Sphinx with:

.. code:: bash

   python -m pip install sphinx

Run the following command in the ``doc/`` folder previously created (where this document is):

.. code:: bash

   sphinx-quickstart

and fill accordingly:

- separate source and build directory: no;
- project name: MadDM;
- authors: list all authors in a “Name Surname” list separated by ``&``;
- project release: leave empty;
- project language: default is ``en``;

Then the folder will be populated by several files.
In particular, we will modify the `conf.py <./conf.py>`__ file as you can check, in order to have dynamic copyright and release, and adding all the various extensions.

Requirements and extensions
---------------------------

Theme
~~~~~

The theme used to build the documentation can be downloaded with

.. code:: bash

   python -m pip install pydata-sphinx-theme

Tutorial and examples gallery
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To automatically build the tutorial and the examples gallery from ipython notebook files containing rst code, we use the `Sphinx-Gallery <https://sphinx-gallery.github.io/stable/index.html>`__ extension:

.. code:: bash

   python -m pip install sphinx-gallery

and we modify the `conf.py <./conf.py>`__ accordingly, by adding the following configuration (or something similar):

.. code:: python

   sphinx_gallery_conf = {
       "examples_dirs": os.path.join(package_root_path, "examples"),
       "gallery_dirs": "tutorial_and_examples",
   }

where:

- ``examples_dirs`` is the path to the folder created in the root of the repository, and containing the ``rst`` notebooks, **this folder should be created**, in this case we will name it ``examples``;
- ``gallery_dirs`` is the name of the folder where the generated files from the extension would be placed after the build.

How to add new examples
~~~~~~~~~~~~~~~~~~~~~~~

The examples are in the folder ``examples`` in the root of the repository.
That folder contains a general header file, which should be named ``GALLERY_HEADER.rst``, that contains some information on the gallery, it will be the landpage for whoever would click on the *Tutorial and examples* entry in the menu.
That file would also show a list of the files present in the gallery.
The gallery is made of python files inside the ``examples`` folder.
See also the `official documentation <https://sphinx-gallery.github.io/stable/getting_started.html#structure-the-examples-folder>`__, but, briefly, there are two kinds of files:

- normal python files, like ``getting_started.py``, they will be parsed and presented in a rich literate programming fashion, but they won’t have any output file/plot associated to them;
- python files with name starting with ``plot_``, they will be executed, their output will be captured and incorporated in the HTML that will be created.

See the `following page <https://sphinx-gallery.github.io/stable/syntax.html#python-script-syntax>`__ for the structure that each of these python files must have.

Updating the documentation
==========================

Whenever a change in the ``.rst`` files is made, the documentation should be rebuilt locally and the ``gh-pages`` branch must be updated accordingly.
The main files to build the documentation are already present when cloning the repository, one just needs to build it with:

.. code:: bash

   make html

And then copy the files created in ``_build`` to the ``gh-pages`` branch.

