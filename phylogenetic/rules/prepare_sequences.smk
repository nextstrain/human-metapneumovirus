"""
This part of the workflow prepares sequences for constructing the phylogenetic tree.

REQUIRED INPUTS:

    metadata    = data/metadata.tsv
    sequences   = data/sequences.fasta
    reference   = data/reference.fasta

OUTPUTS:

    prepared_sequences = results/prepared_sequences.fasta

This part of the workflow usually includes the following steps:

    - augur index
    - augur filter
    - augur align
    - augur mask

See Augur's usage docs for these commands for more details.
"""

rule newreference:
    message:
        """
        Making new reference
        """
    input:
        oldreference = "defaults/reference_{subtype}.gb"
    output:
        newreferencegb = "results/{subtype}/{build}/reference_{build}.gb",
        newreferencefasta = "results/{subtype}/{build}/reference_{build}.fasta",
    params:
        build = lambda w: w.build,
    shell:
        """
        python scripts/newreference.py \
            --reference {input.oldreference} \
            --output-genbank {output.newreferencegb} \
            --output-fasta {output.newreferencefasta} \
            --gene {params.build}
        """

rule filter:
    """
    Filtering to
      - {params.sequences_per_group} sequence(s) per {params.group_by!s}
      - from {params.min_date} onwards
      - excluding strains in {input.exclude}
      - minimum genome length of {params.min_length}
    """
    input:
        sequences = "data/sequences.fasta",
        metadata = "data/metadata.tsv"
    output:
        sequences = "results/{subtype}/{build}/filtered.fasta",
        metadata = "results/{subtype}/{build}/metadata.tsv",
    params:
        group_by = config["filter"]["group_by"],
        sequences_per_group = config["filter"]["sequences_per_group"],
        min_date = config["filter"]["min_date"],
        strain_id = config.get("strain_id_field", "strain"),
        subtype = lambda wc: wc.subtype.upper(),
        build = lambda wc: wc.build,
        min_coverage = lambda w: f'{config["filter"]["min_coverage"].get(w.build, 10000)}',
        subsample_max_sequences = lambda w: config["filter"]["subsample_max_sequences"].get(w.build, 10000),
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --query "subtypes == '{params.subtype}' & {params.build}_coverage > {params.min_coverage}"\
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --group-by {params.group_by} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --min-date {params.min_date}
        """

rule align:
    """
    Aligning sequences to {input.reference}
    """
    input:
        sequences = rules.filter.output.sequences,
        reference = rules.newreference.output.newreferencegb,
    output:
        alignment = "results/{subtype}/{build}/aligned.fasta",
        insertions = "results/{subtype}/{build}/aligned.fasta.insertions.csv",
    threads: 8
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads}
        """

rule filter2:
    """
    Excluding sequences with too many Ns
    """
    input:
        sequences = rules.align.output.alignment,
        metadata = rules.filter.output.metadata,
    output:
        filtered_alignment = "results/{subtype}/{build}/filtered_alignment.fasta",
        filtered_metadata = "results/{subtype}/{build}/filtered_metadata.tsv",
    params:
        min_length = lambda w: config["filter"]["min_length"].get(w.build, 10000),
        strain_id = config.get("strain_id_field", "strain"),
    shell:
        """
        augur filter \
            --metadata-id-columns {params.strain_id} \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --output-sequences {output.filtered_alignment} \
            --output-metadata {output.filtered_metadata}
        """

rule add_insertion:
    """
    Add insertion data to metadata
    """
    input:
        metadata = rules.filter2.output.filtered_metadata,
        insertions = rules.align.output.insertions
    output:
        totalmetadata = 'results/{subtype}/{build}/totalmetadata.tsv'
    params:
        subtype = lambda wc: wc.subtype,
        build = lambda wc: wc.build
    shell:
        """
        python scripts/insertion-metadata.py \
            --aligned-insertions {input.insertions} \
            --metadata {input.metadata}\
            --subtype {params.subtype} \
            --build {params.build} \
            --metadata-insertion {output.totalmetadata}
        """
