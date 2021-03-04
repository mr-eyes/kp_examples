# Create data dir
mkdir data && cd data

# Download and extract the human protein-coding transcripts
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.pc_transcripts.fa.gz
gunzip gencode.v37.pc_transcripts.fa.gz

# Generate names file by gene
grep ">" gencode.v37.pc_transcripts.fa | cut -c2- |  awk -F'|' '{print $0"\t"$2}' > gencode.v37.pc_transcripts.fa.names

