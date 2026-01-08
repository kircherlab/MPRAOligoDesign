.. _Config:

=====================
Config File
=====================

The config file is a yaml file that contains the configuration.You can use one config file per MPRA design. It is divided into :code:`reference` (reference sequence), :code:`datasets` (different input files to design), :code:`tiling` (tiling strategies of predifined regions with variants), and :code:`oligo_design` (filtering and adapters). This is a full example file with all possible configurations. :download:`config/example_config.yaml <../config/example_config.yaml>`.

.. literalinclude:: ../config/example_config.yaml
   :language: yaml
   :linenos:


Note that teh config file is conrolled by json schema. This means that the config file is validated against the schema. If the config file is not valid, the program will exit with an error message. The schema is located in :download:`workflow/schemas/config.schema.yaml <../workflow/schemas/config.schema.yaml>`.

------------------
Reference settings
------------------

The referebce settings are located in the :code:`reference` section. The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_reference
   :end-before: start_datasets

:genome:
    Genome file with lengths of contigs. The full or relative path to the file should be used.
:fasta:
    Genome fasta file (indexed with samtools faidx). The full or relative path to the file should be used.

--------------------
Datasets settings
--------------------

The assignment workflow is configured in the :code:`datasets` section. The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_datasets
   :end-before: start_oligo_length

Each part conntains now a tsv file with a list of samples to include and where the files are located.

:variants_only:
    Samples with only variants (vcf file). Will center the region around the variant.
:regions_only:
    Samples with only regions (bed file). Region length must have the exact size of the oligo design.
:variants_regions:
    Samples with variants and regions (vcf and bed file). Here the tiling strategy takes place (see below)
:sequences_only:
    Only sequences (fasta file). Sequences must have the exact size of the oligo design.

--------------------------------------
Oligo length
--------------------------------------

The experiment workflow is configured in the :code:`oligo_length` section. Each experiment run (contains one experiment file with all replicates of an experiment). The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_oligo_length
   :end-before: end_oligo_length

:end_oligo_length:
    Length of the oligo (excluding adapters).
