rule ancestral:
    input:
        tree=rules.refine.output.tree,
        alignment=rules.exclude.output.filteredsequences, #alignment= rules.align.output.alignment
        annotation=GENBANK_PATH,
        reference='resources/reference.fasta',
    output:
        node_data="output/muts.json",
    params:
        inference = config["ancestral"]["inference"],
        translation_template=r"output/translations/cds_%GENE.translation.fasta",
        output_translation_template=r"output/translations/cds_%GENE.ancestral.fasta",
        genes=" ".join(GENES),
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --inference {params.inference} \
            --annotation {input.annotation} \
            --root-sequence {input.reference} \
            --genes {params.genes} \
            --translations {params.translation_template} \
            --output-node-data {output.node_data} \
            --output-translations {params.output_translation_template}
        """

rule translate:
    """Translating amino acid sequences"""
    input:
        tree = "output/tree.nwk",
        node_data = rules.ancestral.output.node_data,
        reference = "resources/reference.gb",
    output:
        node_data = "output/aa_muts.json"
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
        tree = "output/tree.nwk",
        aa_muts = "output/aa_muts.json",
        nuc_muts = rules.ancestral.output.node_data,
        clades = "resources/clades.tsv",
    output:
        clade_data = "output/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """
