"""
=====================
Use MadDM as a Script
=====================

Run MadDM as a script, allowing for batch processing or automation.

1. Running MadDM as a Script
============================

To run MadDM as a script, you can simply run ``maddm.py`` by providing a text file as an argument containing the
commands you want to execute. For example:

.. code-block:: text

    python bin/maddm.py run.txt

where ``run.txt`` should contain the commands you want to execute, one per line, for example:

.. code-block:: text

    import model DMsimp_s_spin0
    define darkmatter Xr
    generate direct
    output SCAN_example
    launch SCAN_example
    set MXd scan:range(50,700,25)

"""
