.. _Config:

=====================
Config File
=====================

The config file is a yaml file that contains the configuration.You can use one config file per MPRA design. It is divided into :code:`reference` (reference sequence), :code:`datasets` (different input files to design), :code:`tiling` (tiling strategies of predifined regions with variants), and :code:`oligo_design` (filtering and adapters). This is a full example file with all possible configurations. :download:`config/example_config.yml <../config/example_config.yml>`.

.. literalinclude:: ../config/example_config.yml
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
   :end-before: start_tiling

:end_oligo_length:
    Length of the oligo (excluding adapters).


--------------------------------------
Tiling strategy
--------------------------------------

The tiling strategy is configured in the :code:`tiling` section. The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_tiling
   :end-before: start_oligo_design


:remove_edge_variants:
    If set to true, variants that are too close to the edge of an oligo will be removed.
:variant_edge_exclusion:
    Number of bases to exclude variants on edges of oligos when :code:`remove_edge_variants` is true.
:min_overlap:
    Minimum overlap of variant to region to include the variant.
:strategies:
    Different tiling strategies to use. See below for details.

    :centering:
        Center the variant in the oligo. Max defines the maximum length of a region where an oligo should be centered.
    :two_tiles:
        Use two tiles per variant. Max defines the maximum length of a region where two tiles should be used.


--------------------------------------
Oligo design
--------------------------------------

The final oligo design (filtering etc.) is configured in the :code:`oligo_design` section. The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_oligo_design
   :end-before: end_oligo_design


:variants:
    Parameters for variant handling during oligo design

    :use_most_centered_region:
        If variants fits to multiple oligos use the one where it is most centered. Otherwise all are used."
    :remove_unused_regions:    
        Whether to remove regions that are not used for any variant
:filtering:
    Parameters for filtering of oligos

    :max_homopolymer_length:
        Maximum allowed homopolymer length in designed oligos
    :max_simple_repeat_fraction:
        Maximum allowed fraction of simple repeats in designed oligos
:adapters:
    Parameters for adding adapters to oligos
    
    :left:
        Forward adapter sequence to add to the oligos
    :right:
        Reverse adapter sequence to add to the oligos

