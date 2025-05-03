# import some of libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# input the genome data from table browser in UCSC genome browser
hg38_refseqall_genomedata = pd.read_csv("hg38_refseqall_genomedata_20250502.txt", sep="\t")

# sort refseqIDs whose name start with NM_ or NR_
refseqall_nm_or_nr_genomedata = hg38_refseqall_genomedata[hg38_refseqall_genomedata['name'].str.startswith(('NM_', 'NR_'))]
refseqall_nm_genomedata = hg38_refseqall_genomedata[hg38_refseqall_genomedata['name'].str.startswith('NM_')]
refseqall_nr_genomedata = hg38_refseqall_genomedata[hg38_refseqall_genomedata['name'].str.startswith('NR_')]

# sort refseqIDs which have the max of exonCounts
exonCount_idx = refseqall_nm_or_nr_genomedata.groupby("name2")["exonCount"].idxmax()
nm_exonCount_idx = refseqall_nm_genomedata.groupby("name2")["exonCount"].idxmax()
nr_exonCount_idx = refseqall_nr_genomedata.groupby("name2")["exonCount"].idxmax()

refseqall_nm_or_nr_unique_genomedata = refseqall_nm_or_nr_genomedata.loc[exonCount_idx]
refseqall_nm_unique_genomedata = refseqall_nm_genomedata.loc[nm_exonCount_idx]
refseqall_nr_unique_genomedata = refseqall_nr_genomedata.loc[nr_exonCount_idx]

# sort the columns from refseqall_nm_or_nr_unique_genomedata, which we need
# hg38_refseqall_genomedata_columns = ['#bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
# select_features = ['name', 'name2', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds']
select_features = ['name', 'name2', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds']
genome_data = refseqall_nm_or_nr_unique_genomedata[select_features]
nm_genome_data = refseqall_nm_unique_genomedata[select_features]
nr_genome_data = refseqall_nr_unique_genomedata[select_features]

# create intron information
def intron_transform(genome_data):
    # genome_data has exonStarts, exonEnds
    intron_records = []
    
    for _, row in genome_data.iterrows():
        refseq_id = row['name']
        gene = row['name2']
        chrom = row['chrom']
        strand = row['strand']
        exonCount = row['exonCount']
        
        starts = [int(s) for s in row['exonStarts'].strip(',').split(',')] #list内list
        ends = [int(e) for e in row['exonEnds'].strip(',').split(',')] #list内list
        
        for i in range(len(starts) - 1):
            intron_start = ends[i]
            intron_end = starts[i + 1]
            intron_records.append({
                'name': refseq_id,
                'name2': gene,
                'chrom': chrom,
                'strand': strand,
                'exonCount': exonCount,
                'intronCount': exonCount - 1,
                'intron_number': i + 1,
                'intron_start': intron_start,
                'intron_end': intron_end,
                'length': intron_end - intron_start
            })

    # データフレームに変換
    intron_df = pd.DataFrame(intron_records)

    return intron_df

intron_df = intron_transform(genome_data)
intron_nm_df = intron_transform(nm_genome_data)
intron_nr_df = intron_transform(nr_genome_data)

# intron_df.to_csv("intron of all genes including lincRNA and protein coding genes.txt", sep="\t", index=False)
# intron_nm_df.to_csv("intron of protein coding genes.txt", sep="\t", index=False)
# intron_nr_df.to_csv("intron of lincRNA.txt", sep="\t", index=False)

# 描画
# histplot
# x = intron_length, y = counts

def original_hist(intron_df):
    # 平均・中央値を計算
    mean = intron_df["length"].mean()
    median = intron_df["length"].median()
    print(f"{name}のintron長平均値は {mean:.2f}")
    print(f"{name}のintron長中央値は {median:.2f}")
    
    # ヒストグラム描画
    plt.figure(figsize=(8, 5))
    plt.hist(intron_df["length"], bins=50, edgecolor='black') # binの数は増やした方がいいかも
    
    # 平均線（青）と中央値線（赤）を追加
    plt.axvline(mean, color='blue', linestyle='dashed', linewidth=1.0, label=f"Mean: {mean:.0f}")
    plt.axvline(median, color='red', linestyle='dashed', linewidth=1.0, label=f"Median: {median:.0f}")
    
    plt.xlabel("Intron length (bp)")
    plt.ylabel("Count")
    plt.title(f"{name}:Intron Length Distribution")
    plt.grid(True)
    # plt.savefig(f"histgram of {name}")
    plt.show()
    
def cutoff_hist(intron_df, name = 'intron i wanna see'):
    # 外れ値を考慮した（上位5%を除外）
    cutoff = intron_df["length"].quantile(0.95)
    
    # 範囲でフィルタリング
    filtered = intron_df[intron_df["length"] <= cutoff]
    # 平均・中央値を計算
    mean = filtered["length"].mean()
    median = filtered["length"].median()
    print(f"{name}のintron長平均値は {mean:.2f}")
    print(f"{name}のintron長中央値は {median:.2f}")
    
    # ヒストグラム描画
    plt.figure(figsize=(8, 5))
    plt.hist(filtered["length"], bins=50, edgecolor='black')
    # 平均線（青）と中央値線（赤）を追加
    plt.axvline(mean, color='blue', linestyle='dashed', linewidth=1.0, label=f"Mean: {mean:.0f}")
    plt.axvline(median, color='red', linestyle='dashed', linewidth=1.0, label=f"Median: {median:.0f}")
    
    plt.xlabel("Intron length (bp)")
    plt.ylabel("Count")
    plt.title(f"{name}:Intron Length Distribution (≤ {int(cutoff)} bp)")
    plt.grid(True)
    # plt.savefig(f"histgram of {name}")
    plt.show()
