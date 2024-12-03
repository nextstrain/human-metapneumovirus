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
        ref = "defaults/reference.fasta",
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
# rule subtypes:
#     input:
#         metadata = rules.extend_metadata.output.metadata,
#         sequences = "data/sequences.fasta",
#     output:
#         metadata = "data/metadata_{subtype}",
#         sequences = "data/sequences_{subtype}",
#     params:
#         subtype = lambda wc: wc.subtype
#     shell:
#         """
#         augur filter \
#             --metadata {input.metadata} \
#             --sequences {input.sequences} \
#             --query "subtypes == '{params.subtype}'" \
#             --output-sequences {output.sequences} \
#             --output-metadata {output.metadata}
#         """