# intron_informatics_2025

This project analyzes intron information based on RefSeq data.

## Data Source

The raw data was downloaded from the [UCSC Genome Browser](https://genome.ucsc.edu/) with the following settings:

- **Genome**: Human (Homo sapiens)  
- **Assembly**: Dec. 2013 (GRCh38/hg38)  
- **Group**: Genes and Gene Prediction  
- **Track**: NCBI RefSeq  
- **Table**: RefSeq All  
- **Region**: Whole Genome

## histgram of intron length
外れ値上位5%を削除しています

### all genes (NM_(protein) and NR_(lincRNA)
all intronのintron長平均値（青）は 3072.61  
all intronのintron長中央値（赤）は 1430.00  
<img src="/images/histgram of all intron.png" width="500">

### protein coding genes, NM_(protein)
intron of genes which code proteinのintron長平均値（青）は 2886.95  
intron of genes which code proteinのintron長中央値（赤）は 1395.50  
<img src="/images/histgram of intron of genes which code protein.png" width="500">

intronlengthを0から2000の間に絞った図  
<img src="/images/histgram_sperange of intron of genes which code protein range=0-2000.png" width="500">

### lincRNA genes, NR_(lincRNA)
intron of genes which code lincRNAのintron長平均値（青）は 3868.43  
intron of genes which code lincRNAのintron長中央値（赤）は 1651.50  
<img src="/images/histgram of intron of genes which code lincRNA.png" width="500">

### protein coding genes, NM_(protein).ver2 #未定義な配列を除去 #上位3%のlengthを持つレコードを除去 #イントロンの個数が最も多いようなrefseqIDを取得
intron of genes which code proteinのintron長平均値（青）は 3423.61  
intron of genes which code proteinのintron長中央値（赤）は 1454.00  
<img src="/images/histgram_20250507_protein coding genes.png" width="500">

intronlengthを0から2000の間に絞った図  
<img src="/images/histgram_20250507_protein coding genes_0-2000.png" width="500">
