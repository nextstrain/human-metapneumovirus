### add insertion information to metadata ###
import pandas as pd
import csv

df1 = pd.read_csv ("~/mpv/mpv-nextstrain-workflow/phylogenetic/data/metadata.tsv", sep = '\t')
df2 = pd.read_csv("~/mpv/mpv-nextstrain-workflow/phylogenetic/results/metadata-insertion-6738.csv")

# Perform a left merge
totalmetadata = pd.merge(df1, df2, on='accession', how='left')
totalmetadata['boolean'] = totalmetadata['boolean'].fillna(False)
totalmetadata.rename(columns={'boolean': 'insertion6738'}, inplace=True)
totalmetadata.to_csv("~/mpv/mpv-nextstrain-workflow/phylogenetic/results/totalmetadata.tsv", sep='\t', index= False)