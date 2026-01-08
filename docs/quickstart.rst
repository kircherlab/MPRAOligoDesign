.. _Getting started:


=====================
Getting started
=====================


To run MPRAOligoDesign, first activate the snakemake environment, e-g., via:

.. code-block:: bash

    conda activate snakemake

All inputs of the workflow are configured via a config file. For more information on the parameters see :doc:`config`.



And then run the main workflow with:

.. code-block:: bash

    snakemake --use-conda --cores $N --configfile config/config.yaml
