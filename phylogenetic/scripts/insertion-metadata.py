import pandas as pd
import csv
import re

def add_insertions(subtype, build, metadata, alignedinsertions, metadata_insertion):

    alignedinsertionsdf = pd.read_csv(alignedinsertions, index_col="strain")
    metadatadf = pd.read_csv(metadata, sep = '\t', index_col="accession")

    if subtype == "A" and build != "F":

        if build =="genome":
            min_position = 6247
            max_position = 6956

        if build == "G":
            min_position = 1
            max_position = 710 

        #select columns in gene G with insertions
        column_with_insertion = []
        for column in alignedinsertionsdf.columns:
            match = re.match(r"^insertion.*pos (\d+)", column)
            if match is None:
                continue
            if min_position <= int(match.group(1)) <= max_position:
                column_with_insertion.append(column)

        #select strains in gene G with insertions longer than 30 and add column with insertion data to metadata file
        index_with_insertion = []
        for index, row in alignedinsertionsdf.iterrows():
            has_insertion = any(len(str(row[col])) > 30 for col in column_with_insertion)
            if has_insertion:
                index_with_insertion.append(index)
            print(index_with_insertion)
        metadatadf['insertion'] = metadatadf.index.isin(index_with_insertion)
        metadatadf.to_csv(metadata_insertion, sep='\t', index= "accession")

    elif subtype == "B" or build == "F":
        metadatadf.to_csv(metadata_insertion, sep='\t', index= "accession")

if __name__ == '__main__':
    import argparse

    parser = parser = argparse.ArgumentParser(description='parse metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--aligned_insertions', help="input aligned insertions file")           
    parser.add_argument('--metadata', help="input metadata file")
    parser.add_argument('--metadata_insertion', help="final output metadata with insertion in G gene")
    parser.add_argument('--subtype')
    parser.add_argument('--build')
    args = parser.parse_args()

    add_insertions(args.subtype, args.build, args.metadata, args.aligned_insertions, args.metadata_insertion)

