import pandas as pd
import csv

if __name__ == '__main__':
    import argparse

    parser = parser = argparse.ArgumentParser(description='parse metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--alignedinsertions', help="input aligned insertions file")           
    parser.add_argument('--metadata', help="input metadata file")
    parser.add_argument('--metadata_insertion', help="final output metadata with >30 bp insertion in G gene")
    args = parser.parse_args()

    ##################################################
    #   Part 1
    ##################################################
    
    alignedinsertionsdf = pd.read_csv(args.alignedinsertions, index_col=False)
    metadatadf = pd.read_csv(args.metadata, sep = '\t', index_col=False)

    #create a new file with accession number and insertion metadata
    d = {'strain': alignedinsertionsdf['strain'].values, 'boolean': alignedinsertionsdf['insertion: 180bp @ ref pos 6738'].str.len() > 30}
    df = pd.DataFrame(data=d)
    df.rename(columns={'strain': 'accession'}, inplace=True)

    #merge new file with insertion data and metadata file
    totalmetadata = pd.merge(metadatadf, df, on='accession', how='left')
    totalmetadata['boolean'] = totalmetadata['boolean'].fillna(False)
    totalmetadata.rename(columns={'boolean': 'insertion6738'}, inplace=True)
    totalmetadata.to_csv(args.metadata_insertion, sep='\t', index= False)
