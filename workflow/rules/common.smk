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


def isVariantsAndRegionsSample(sample):
    """
    Returns True if the sample is a variants and regions sample, False otherwise.
    """
    return "bed_file" in datasets.columns and not pd.isnull(datasets.loc[sample, "bed_file"])

def getRegionDatasets():
    """
    Returns a list of samples that contain regions.
    """
    if "bed_file" in datasets.columns:
        return datasets[not datasets["bed_file"].isnull()].index.tolist()
    else:
        return[]

def getVariantDatasets():
    """
    Returns a list of samples that contain variants.
    """
    if "vcf_file" in datasets.columns:
        return datasets[not datasets["vcf_file"].isnull()].index.tolist()
    else:
        []
