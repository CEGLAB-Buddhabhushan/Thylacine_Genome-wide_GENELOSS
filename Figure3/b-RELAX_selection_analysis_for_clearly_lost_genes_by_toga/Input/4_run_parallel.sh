cat transcript.lst | parallel -j 16 bash 3_extract_seq_chain_toga.sh {}
