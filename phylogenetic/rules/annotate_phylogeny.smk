"""
This part of the workflow creates additonal annotations for the phylogenetic tree.

REQUIRED INPUTS:

    metadata            = data/metadata.tsv
    prepared_sequences  = results/prepared_sequences.fasta
    tree                = results/tree.nwk

OUTPUTS:

    node_data = results/*.json

    There are no required outputs for this part of the workflow as it depends
    on which annotations are created. All outputs are expected to be node data
    JSON files that can be fed into `augur export`.

    See Nextstrain's data format docs for more details on node data JSONs:
    https://docs.nextstrain.org/page/reference/data-formats.html

This part of the workflow usually includes the following steps:

    - augur traits
    - augur ancestral
    - augur translate
    - augur clades

See Augur's usage docs for these commands for more details.

Custom node data files can also be produced by build-specific scripts in addition
to the ones produced by Augur commands.
"""
rule ancestral:
    """Reconstructing ancestral sequences and mutations"""
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment,
    output:
        node_data = "results/nt_muts_{subtype}.json"
    params:
        inference = config["ancestral"]["inference"]
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    """Translating amino acid sequences"""
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = "defaults/reference_{subtype}.gb",
    output:
        node_data = "results/aa_muts_{subtype}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule clades:
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = "defaults/clades_{subtype}.tsv",
    output:
        clade_data = "results/clades_{subtype}.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

rule traits:
    """
    Inferring ancestral traits for {params.columns!s}
      - increase uncertainty of reconstruction by {params.sampling_bias_correction} to partially account for sampling bias
    """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.add_insertion.output.totalmetadata,
    output:
        node_data = "results/traits_{subtype}.json",
    params:
        columns = config["traits"]["columns"],
        sampling_bias_correction = config["traits"]["sampling_bias_correction"],
        strain_id = config.get("strain_id_field", "strain"),
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction}
        """
