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
    define darkmatter xd
    generate direct
    output SCAN_example
    launch SCAN_example
    direct_nucleon = ON
    set limit
    set direct_nucleon_targets xe, ar
    set MXd scan:range(50,700,25)

After ``launch``, the lines answer the launch prompt: switches (``direct_nucleon = ON``, ``relic = ON``, ...), ``set`` commands
for the cards and shortcuts such as ``set limit`` or ``set limit_and_rates`` for the direct detection mode. The run starts at the
end of the file, as if Enter was pressed. Note that the switch names are those shown at the launch prompt: ``direct_nucleon`` and
``direct_electron`` for direct detection (``direct`` is only the name of the process in ``generate direct``).

"""

# sphinx_gallery_thumbnail_path = '_static/thumbnail/scripts_thumb.png'