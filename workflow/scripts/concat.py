import pandas as pd

frames = list(snakemake.input)
dfs = []
for f in frames:
    dfs.append(pd.read_csv(f, sep='\t'))

pd.concat(dfs).to_csv(str(snakemake.output), sep='\t', index=False)
