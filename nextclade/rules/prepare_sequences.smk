REFERENCE_ACCESSION = "NC_039199.1"
TAXON_ID = 162145
GENES = ["N", "P", "M", "F", "M2-1", "M2-2", "SH", "G", "L"]
ALLOWED_DIVERGENCE = "4000"
GFF_PATH = "resources/genome_annotation.ggf3"
GENBANK_PATH = "resources/reference.gb"
REFERENCE_PATH = "resources/reference.fasta"
README_PATH = "resources/README.md"
CHANGELOG_PATH = "resources/CHANGELOG.md"
ROOTING = "mid_point"  # alternative root using outgroup, e.g. the reference NC_039199.1
MIN_LENGTH = 1000 # Minimal length of sequences for example dataset.


# #rule add_reference_to_include:
#     """
#     Create an include file for augur filter
#     """
#     #input:
#     #    include="resources/include.txt",
#     #output:
#     #    "results/include.txt",
#     #shell:
#         """
#         echo "{REFERENCE_ACCESSION}" >> results/include.txt
#         """

rule filter:
    """
    Filtering to
      - {params.sequences_per_group} sequence(s) per {params.group_by!s}
      - from {params.min_date} onwards
      - excluding strains in {input.exclude}
      - minimum genome length of {params.min_length} 
    """
    input:
        sequences = "resources/sequences.fasta",
        metadata = "resources/metadata.tsv",
    output:
        sequences = "output/filtered.fasta",
    params:
        group_by = config["filter"]["group_by"],
        sequences_per_group = config["filter"]["sequences_per_group"],
        min_date = config["filter"]["min_date"],
        min_length = config["filter"]["min_length"],
        strain_id = config.get("strain_id_field", "strain"),
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length}
        """

rule align:
    """
    Aligning sequences to {input.reference}
    """
    input:
        sequences = rules.filter.output.sequences,
        reference = "resources/reference.fasta",
        pathogen_json = "resources/pathogen.json",
        annotations = GFF_PATH,
    output:
        alignment = "output/aligned.fasta",
        tsv = "output/nextclade.tsv",
    params:
        translation_template=lambda w: "output/translations/cds_{cds}.translation.fasta",
    shell:
        """
        nextclade3 run \
            {input.sequences} \
            --input-ref {input.reference} \
            --input-pathogen-json {input.pathogen_json} \
            --input-annotation {input.annotations} \
            --include-reference \
            --output-translations {params.translation_template} \
            --output-tsv {output.tsv} \
            --output-fasta {output.alignment}
        """


rule get_outliers:
    """
    Automatically identify sequences with >{ALLOWED_DIVERGENCE} substitutions
    (likely to be sequencing errors or low quality/misannotated sequences) and put them in outliers.txt
    """
    input:
        nextclade= rules.align.output.tsv,
    output:
        outliers="output/outliers.txt",
        tmp="tmp/outliers.txt",
    params:
        allowed_divergence=lambda w: ALLOWED_DIVERGENCE,
    shell:
        """
        tsv-filter -H -v --is-numeric totalSubstitutions {input.nextclade} \
        > {output.tmp}
        tsv-filter -H \
            --is-numeric totalSubstitutions \
            --gt totalSubstitutions:{params.allowed_divergence} \
            {input.nextclade} \
        | tail -n +2 >> {output.tmp}
        cat {output.tmp} \
        | tsv-select -H -f seqName \
        | tail -n +2 > {output.outliers}
        """

rule short_sequences:
    """
    Identify sequences shorter than {MIN_LENGTH} and put them in short_sequences.txt
    """
    input:
        sequences = "resources/sequences.fasta",
    params:
        min_length = MIN_LENGTH
    output:
        short_sequences = "output/short_sequences.txt"
    shell:
        """
        seqkit seq {input.sequences} -M {params.min_length} -n > {output.short_sequences}
        """

rule exclude:
    """
    Rule to allow for manual and automatic exclusion of sequences
    without triggering a new subsampling that could
    surface new bad sequences resulting in an infinite loop
    """
    input:
        sequences=rules.align.output.alignment,
        metadata="resources/metadata.tsv",
        exclude="resources/exclude.txt",
        outliers="output/outliers.txt",
    params:
        strain_id = config.get("strain_id_field", "strain"),
    output:
        filteredsequences="output/filteredaligned.fasta",
        filteredmetadata="output/filteredmetadata.tsv",
        strains="output/treestrains.txt",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --exclude {input.exclude} {input.outliers} \
            --output {output.filteredsequences} \
            --output-metadata {output.filteredmetadata} \
            --output-strains {output.strains}
        """