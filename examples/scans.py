"""

==========================
Scan over model parameters
==========================

Scan over a specified range of model parameters.

1. Running a Scan
=================

First, you need to set up your model and process as usual. For example, if you want to scan over a simplified scalar dark matter model, you can do the following:

.. code-block:: text

    MadDM> import model DMsimp_s_spin0
    MadDM> define darkmatter Xr
    MadDM> generate direct
    MadDM> output SCAN_example
    MadDM> launch SCAN_example

At the launch prompt, set the parameters you want to scan, using ``set MXd scan:`` followed by a python array or ``range`` function.
For instance, if you want to scan the dark matter mass (MXd) from 50 GeV to <700 GeV in steps of 25 GeV, you can do:

.. code-block:: text

    MadDM> set MXd scan:range(50,700,25)

If you want to scan only over some specific values, you can also use python arrays, as follows:

.. code-block:: text

    MadDM> set MXd scan:[10,20,40,80,160,320,640]

or you can use list comprehension. For example, to scan over a logarithmic scale from 1 GeV to 10 GeV you can do:

.. code-block:: text

    MadDM> set MXd scan:[10 ** (i * 0.01) for i in range(101)]

Then press Enter to start the scan. For scans, you will find in the ``output`` directory a ``run_XX_YY`` folder containing the ``maddm.out`` file  for each ``YY`` iteration.
You will also find a ``scan_run_XX.txt`` file that summarizes the scan parameters and results.
You can also make a scan over multiple parameters at once. For example, if you want to scan over the dark matter mass ``MXd`` and the coupling ,``gsxd`` you can do:

.. code-block:: text

    MadDM> set MXd scan:range(50,700,25)
    MadDM> set gsxd scan:range(10,100,10)

Then press Enter to start the scan. MadDM will generate one run for each combination of the parameters you specified, iterating first over the last parameter in the `param_card.dat`.
Once all values of that parameter are exhausted, it steps the second-to-last parameter and repeats the process, and so on, like nested ``for`` loops, starting from the outermost (first)
parameter to the innermost (last).

Exiting MadDM
=============

.. code-block:: text

    MadDM> quit

"""

# sphinx_gallery_thumbnail_path = '_static/thumbnail/scans_thumb.png'