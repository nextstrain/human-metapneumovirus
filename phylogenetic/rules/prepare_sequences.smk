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
        metadata = "data/metadata.tsv",
    output:
        sequences = "results/filtered.fasta"
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
        sequences = "results/filtered.fasta",
        reference = "defaults/reference.gb",
    output:
        alignment = "results/aligned.fasta",
        insertions = "results/aligned.fasta.insertions.csv",
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
        metadata = "data/metadata.tsv",
        insertions = rules.align.output.insertions
    output: 
        totalmetadata = 'results/totalmetadata.tsv'
    shell:
        """
        python scripts/insertion-metadata.py \
            --alignedinsertions {input.insertions} \
            --metadata {input.metadata}\
            --metadata_insertion {output.totalmetadata}
        """