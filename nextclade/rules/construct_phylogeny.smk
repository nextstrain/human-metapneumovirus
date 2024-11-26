rule tree:
    input:
        alignment=rules.exclude.output.filteredsequences, 
    output:
        tree="output/tree_raw.nwk",
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
        """

rule root:
    input:
        tree= rules.tree.output.tree
    output:
        tree = "output/tree_rooted.nwk"
    run:
        from Bio import Phylo
        T = Phylo.read(input.tree, 'newick')
        T.root_at_midpoint()
        Phylo.write(T, output.tree, 'newick')

rule refine:
    input:
        tree= rules.root.output.tree,
        alignment=rules.exclude.output.filteredsequences, 
    output:
        tree="output/tree.nwk",
        node_data="output/branch_lengths.json",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --keep-root \
            --keep-polytomies \
            --output-node-data {output.node_data} \
            --output-tree {output.tree}
        """
