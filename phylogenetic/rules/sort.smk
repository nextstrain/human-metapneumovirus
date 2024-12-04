"""
This part of the workflow handles sorting downloaded sequences and metadata
into a and b by aligning them to reference sequences.

It produces output files as

    metadata = "data/{type}/metadata.tsv"

    sequences = "data/{type}/sequences.fasta"

"""

rule nextclade:
    input:
        sequences = "data/sequences.fasta",
        ref = "defaults/reference_A.fasta",
        dataset = "../nextclade/dataset"
    output:
        nextclade = "data/nextclade.tsv"
    params:
        output_columns = "seqName clade qc.overallScore qc.overallStatus alignmentScore  alignmentStart  alignmentEnd  coverage dynamic"
    threads: 8
    shell:
        """
        nextclade3 run -D {input.dataset}  -j {threads} \
                          --output-columns-selection {params.output_columns} \
                          --output-tsv {output.nextclade} \
                          {input.sequences}
        """

rule extend_metadata:
    input:
        nextclade = rules.nextclade.output.nextclade,
        metadata = "data/metadata.tsv"
    output:
        metadata = "data/extended_metadata.tsv"
    shell:
        """
        python3 scripts/extend-metadata.py --metadata {input.metadata} \
                                       --id-field accession \
                                       --nextclade {input.nextclade} \
                                       --output {output.metadata}
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
