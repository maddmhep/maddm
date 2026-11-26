"""

==========================
Scan over model parameters
==========================

Scan over a specified range of model parameters.

1. Running a Scan
=================

First, you need to set up your model and process as usual. For example, to scan the DM-electron scattering of a simplified scalar-mediated
dark matter model, you can do the following:

.. code-block:: text

    MadDM> import model eDMsimp_s_spin0
    MadDM> define darkmatter xd
    MadDM> generate direct
    MadDM> output SCAN_example
    MadDM> launch SCAN_example

This model has only DM-electron couplings, so at the launch prompt only ``direct_electron`` is ON (``direct_nucleon`` is OFF).
At the launch prompt, set the parameters you want to scan, using ``set MXd scan:`` followed by a python array or ``range`` function.
For instance, if you want to scan the dark matter mass (MXd) from 50 GeV to <700 GeV in steps of 25 GeV, you can do:

.. code-block:: text

    MadDM> set MXd scan:range(50,700,25)

Note that ``direct_electron`` skips DM masses above 1 GeV unless ``direct_electron_mode`` is set to ``always``
(``set direct_electron_mode always``): for DM-electron scattering a sub-GeV range, as in the examples below, is the relevant one.

If you want to scan only over some specific values, you can also use python arrays, as follows:

.. code-block:: text

    MadDM> set MXd scan:[0.01,0.02,0.04,0.10,0.20,0.4,0.8,1]

or you can use list comprehension. For example, to scan over a logarithmic scale from 10 MeV to 1 GeV you can do:

.. code-block:: text

    MadDM> set MXd scan:[10 ** (-2 + i * 0.2) for i in range(11)]

If you only need the cross sections and the exclusion limits, and not the recoil spectra of each point, set ``direct_mode`` to ``limit``
before starting the scan, which is faster (``set limit_and_rates``, the default, also computes the spectra and, for DM-electron
scattering, the XENON10 and XENON1T results):

.. code-block:: text

    MadDM> set limit

Then press Enter to start the scan. For scans, you will find in the ``output`` directory a ``run_XX_YY`` folder for each ``YY`` iteration,
with the ``MadDM_results.txt`` and ``maddm.out`` files and the output files kept by ``save_output``.
You will also find a ``scan_run_XX.txt`` file that summarizes the scan parameters and results.
You can also make a scan over multiple parameters at once. For example, if you want to scan over the dark matter mass ``MXd`` and the coupling ``gsxd``, you can do:

.. code-block:: text

    MadDM> set MXd scan:[10 ** (-2 + i * 0.2) for i in range(11)]
    MadDM> set gsxd scan:range(100,1000,100)

Then press Enter to start the scan. MadDM will generate one run for each combination of the parameters you specified, iterating first over the last parameter in the `param_card.dat`.
Once all values of that parameter are exhausted, it steps the second-to-last parameter and repeats the process, and so on, like nested ``for`` loops, starting from the outermost (first)
parameter to the innermost (last).

You will find the results of the scan in ``SCAN_example/output/scan_run_01.txt``. In this example, this file will contain,
for each value of the dark matter mass and coupling you scanned, the reference DM-electron cross section ``sigma_e``, the XENONnT
limits ``lim_sigma_e`` and ``lim_sigma_e_rescaled``, the DAMIC-M limits ``lim_sigma_e_DAMICM`` (published) and
``lim_sigma_e_DAMICM_rescaled`` (rescaled to the RAPIDD silicon response), the DarkSide-50 limit ``lim_sigma_e_DS50`` (-1 outside
their mass range), with ``limit_and_rates`` the XENON10 and XENON1T limits and p-values ``lim_sigma_e_Xenon10``, ``pvalue_Xenon10``,
``lim_sigma_e_Xenon1T``, ``pvalue_Xenon1T`` (the only tests computed from the recoil spectrum of the model), and the
heavy-mediator flag ``dd_e_light_mediator_warning`` (0 if the contact approximation holds at that point). For DM-nucleon scattering (``direct_nucleon``) the file contains the cross sections ``sigmaN_*``, the
equivalent cross sections ``sigmaN_*_eq``, their limits ``lim_sigmaN_*`` and the flag ``dd_light_mediator_warning``.

Exiting MadDM
=============

.. code-block:: text

    MadDM> quit

"""

# sphinx_gallery_thumbnail_path = '_static/thumbnail/scans_thumb.png'