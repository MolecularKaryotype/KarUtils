import pandas as pd


df = pd.read_csv('acrocentric_telo_cen.bed', sep='\t')
df['len'] = df['EndPos'] - df['StartPos']
df = df[df['Type'].isin(['acrocentric-telomere1'])]
print(df['len'].sum())