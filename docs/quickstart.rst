.. _Getting started:


=====================
Getting started
=====================


To run MPRAOligoDesign, first activate the snakemake environment, e-g., via:

.. code-block:: bash

    conda activate snakemake

All inputs of the workflow are configured via a config file. For more information on the parameters see :doc:`config`.

Pre required files
--------------------

Before running the workflow we need the human reference genome from UCSC (GRCh38 build) and the genome size file. Be aware this is around 3GB large after extraction. You can download them with:

.. code-block:: bash

    wget -P resources/reference https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
    wget -P resources/reference https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.fai
    wget -P resources/reference https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes

    gunzip resources/reference/hg38.fa.gz

You need to set the correct path in the config file. Here you have to modify the entries ``genome`` and ``fasta`` under ``reference``. We will use the config file ``config/config_quickstart.yml`` in this example.

Running the workflow
------------------------

We can combine 4 different designing routines. ``variants_only``, ``regions_only``, ``variants_regions``, and ``sequences_only``. See :doc:`config` for more information on them. In this quickstart we will use all four routines. The input files are already provided in the ``data/quickstart`` folder. You can run the workflow with:

.. code-block:: bash

    snakemake --cores 1 --configfile config/config_quickstart.yml

This should design oligogs for all 3 routines. The results will be stored in the ``results/final_design`` folder. 
