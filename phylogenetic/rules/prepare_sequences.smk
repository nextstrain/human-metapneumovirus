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
        metadata = rules.extend_metadata.output.metadata,
    output:
        sequences = "results/filtered_{subtype}.fasta",
        metadata = "results/metadata_{subtype}.tsv",
    params:
        group_by = config["filter"]["group_by"],
        sequences_per_group = config["filter"]["sequences_per_group"],
        min_date = config["filter"]["min_date"],
        min_length = config["filter"]["min_length"],
        strain_id = config.get("strain_id_field", "strain"),
        subtype = lambda wc: wc.subtype
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --query "subtypes == '{params.subtype}'" \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
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
        reference = "defaults/reference_{subtype}.gb",
    output:
        alignment = "results/aligned_{subtype}.fasta",
        insertions = "results/aligned_{subtype}.fasta.insertions.csv",
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference
        """
####################################
# add insertion data to metadata
####################################

rule add_insertion:
    input:
        metadata = rules.filter.output.metadata,
        insertions = rules.align.output.insertions
    output: 
        totalmetadata = 'results/totalmetadata_{subtype}.tsv'
    params:
        subtype = lambda wc: wc.subtype
    shell:
        """
        python scripts/insertion-metadata.py \
            --alignedinsertions {input.insertions} \
            --metadata {input.metadata}\
            --subtype {params.subtype} \
            --metadata_insertion {output.totalmetadata}
        """