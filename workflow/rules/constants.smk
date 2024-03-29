HAPS = ["hap1", "hap2"]
TECH = config.get("tech", "ont")

# --------  Constants: PARAMS -------- #
MINIMAP2_PARAMS = config["minimap2"].get("params", "")

# --------  Constants: VERSION -------- #
SAMTOOLS_VERSION = config["samtools"]["version"]
BCFTOOLS_VERSION = config["bcftools"]["version"]
SNIFFLES_VERSION = config["sniffles"]["version"]
LONGPHASE_VERSION = config["longphase"]["version"]
SAMBAMBA_VERSION = config["sambamba"]["version"]
MINIMAP2_VERSION = config["minimap2"]["version"]
METHYLINK_VERSION = config["methylink"]["version"]
SEQTK_VERSION = config["seqtk"]["version"]
CANU_VERSION = config["canu"]["version"]
MODKIT_VERSION = config["modkit"]["version"]
BEDTOOLS_VERSION = config["bedtools"]["version"]
YAK_VERSION = config["yak"]["version"]
HIFIASM_VERSION = config["hifiasm"]["version"]
pbCpGtools_VERSION = config["pb-CpG-tools"]["version"]
UCSC_VERSION = config["ucsc"]["version"]

# --------  Constants: CONTAINERS -------- #
CLAIR3_CNTR = config["clair3"]["container"]
