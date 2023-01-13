################################
#### Global functions       ####
################################
from snakemake.workflow import srcdir

SCRIPTS_DIR = srcdir("../scripts")


def getScript(name):
    return "%s/%s" % (SCRIPTS_DIR, name)


REFERENCE_DIR = srcdir("../../reference")


def getReference(name):
    return "%s/%s" % (REFERENCE_DIR, name)


from snakemake.utils import validate
import pandas as pd


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


##### load config and sample sheets #####

# preferrred to use --configfile instead of hard-coded config file
# configfile: "config/config.yaml"

validate(config, schema="../schemas/config.schema.yaml")

datasets = pd.DataFrame()

if "variants_regions" in config["datasets"]:
    df = pd.read_csv(config["datasets"]["variants_regions"], sep="\t").set_index(
        "sample", drop=False
    )
    df.index.names = ["sample_id"]
    validate(df, schema="../schemas/datasets_variants_regions.schema.yaml")
    datasets = pd.concat([datasets, df])
if "variants_only" in config["datasets"]:
    df = pd.read_csv(config["datasets"]["variants_only"], sep="\t").set_index(
        "sample", drop=False
    )
    df.index.names = ["sample_id"]
    validate(df, schema="../schemas/datasets_variants_only.schema.yaml")
    datasets = pd.concat([datasets, df])
if "regions_only" in config["datasets"]:
    df = pd.read_csv(config["datasets"]["regions_only"], sep="\t").set_index(
        "sample", drop=False
    )
    df.index.names = ["sample_id"]
    validate(df, schema="../schemas/datasets_regions_only.schema.yaml")
    datasets = pd.concat([datasets, df])
if "sequences_only" in config["datasets"]:
    df = pd.read_csv(config["datasets"]["sequences_only"], sep="\t").set_index(
        "sample", drop=False
    )
    df.index.names = ["sample_id"]
    validate(df, schema="../schemas/datasets_sequences_only.schema.yaml")
    datasets = pd.concat([datasets, df])


def isVariantsAndRegionsSample(sample):
    """
    Returns True if the sample is a variants and regions sample, False otherwise.
    """
    return (
        "vcf_file" in datasets.columns and pd.notnull(datasets.loc[sample, "vcf_file"])
    ) and (
        "bed_file" in datasets.columns and pd.notnull(datasets.loc[sample, "bed_file"])
    )


def isVariantsOnlySample(sample):
    """
    Returns True if the sample is a variants and regions sample, False otherwise.
    """
    return (
        "vcf_file" in datasets.columns and pd.notnull(datasets.loc[sample, "vcf_file"])
    ) and not isVariantsAndRegionsSample(sample)


def isRegionsOnlySample(sample):
    """
    Returns True if the sample is a regions sample, False otherwise.
    """
    return (
        "bed_file" in datasets.columns
        and pd.notnull(datasets.loc[sample, "bed_file"])
        and not isVariantsAndRegionsSample(sample)
    )


def isSequenceOnlySample(sample):
    """
    Returns True if the sample is a sequence only sample, False otherwise.
    """

    return "fasta_file" in datasets.columns and pd.notnull(
        datasets.loc[sample, "fasta_file"]
    )


def getRegionDatasets():
    """
    Returns a list of samples that contain regions.
    """
    if "bed_file" in datasets.columns:
        return datasets[datasets["bed_file"].notnull()].index.tolist()
    else:
        return []


def getVariantDatasets():
    """
    Returns a list of samples that contain variants.
    """
    if "vcf_file" in datasets.columns:
        return datasets[datasets["vcf_file"].notnull()].index.tolist()
    else:
        return []


def getFinalRegionFile(sample):
    if isVariantsAndRegionsSample(sample) or isVariantsOnlySample(sample):
        print("R Variants", sample)
        return "results/oligo_design/{sample}/design_variants_filtered.regions.bed.gz"
    elif isRegionsOnlySample(sample):
        print("R RegionsOnly", sample)
        return "results/oligo_design/{sample}/design_regions_filtered.regions.bed.gz"


def getFinalDesignFile(sample):
    if isVariantsAndRegionsSample(sample) or isVariantsOnlySample(sample):
        print("D Variants", sample)
        return "results/oligo_design/{sample}/design_variants_filtered.design.fa"
    elif isRegionsOnlySample(sample):
        print("D RegionsOnly", sample)
        return "results/oligo_design/{sample}/design_regions_filtered.design.fa"
    elif isSequenceOnlySample(sample):
        print("D SequenceOnly", sample)
        return ("results/oligo_design/{sample}/design_sequences_filtered.design.fa",)
