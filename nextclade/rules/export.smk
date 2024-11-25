rule export:
    input:
        tree="output/tree.nwk",
        metadata="output/filteredmetadata.tsv",
        mutations="output/muts.json",
        branch_lengths="output/branch_lengths.json",
        clades="output/clades.json",
        auspice_config="resources/auspice_config.json",
    output:
        auspice="output/auspice.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --auspice-config {input.auspice_config} \
            --node-data {input.mutations} {input.branch_lengths} {input.clades} \
            --output {output.auspice}
        """


rule subsample_example_sequences:
    input:
        all_sequences="resources/sequences.fasta",
        tree_strains="output/treestrains.txt",
        short_sequences = rules.short_sequences.output.short_sequences
    output:
        exclude ="output/exclude.txt",
        example_sequences="output/example_sequences.fasta",
    shell:
       """
        # Exclude short and outlier sequences from all sequences
        cat {input.tree_strains} {input.short_sequences} > {output.exclude} 
        seqkit grep -v -f {output.exclude} {input.all_sequences} \
        | seqkit grep -v -f {input.short_sequences} \
        | seqkit sample -n 100 -s 42 > {output.example_sequences}
        """

rule assemble_dataset:
    input:
        tree=rules.export.output.auspice,
        reference=REFERENCE_PATH,
        annotation=GFF_PATH,
        sequences=rules.subsample_example_sequences.output.example_sequences,
        pathogen="resources/pathogen.json",
        readme= README_PATH,
        changelog=CHANGELOG_PATH,
    output:
        tree="dataset/tree.json",
        reference="dataset/reference.fasta",
        annotation="dataset/genome_annotation.gff3",
        sequences="dataset/sequences.fasta",
        pathogen="dataset/pathogen.json",
        readme="dataset/README.md",
        changelog="dataset/CHANGELOG.md",
        dataset_zip="dataset.zip",
    shell:
        """
        cp {input.tree} {output.tree}
        cp {input.reference} {output.reference}
        cp {input.annotation} {output.annotation}
        cp {input.sequences} {output.sequences}
        cp {input.pathogen} {output.pathogen}
        cp {input.readme} {output.readme}
        cp {input.changelog} {output.changelog}
        zip -rj dataset.zip  dataset/*
        """


rule test:
    input:
        dataset="dataset.zip",
        sequences="dataset/sequences.fasta",
    output:
        output=directory("test_out"),
    shell:
        """
        nextclade3 run \
            --input-dataset {input.dataset} \
            --output-all {output.output} \
            {input.sequences}
        """