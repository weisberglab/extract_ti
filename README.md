# Extract Ti/Ri plasmid and genes
Extract Ti/Ri oncogenic plasmid based on BEAV output

# What is does
- Look for all known Ti/Ri genes in the GenBank file generated by BEAV
- Extract individual Ti/Ri genes and save them in `.ffa` & `.faa` format
- Summarize known Ti/Ri genes predicted by BEAV in  `.gff3` format
- Extract the whole oncogenic plasmid(s) in a seperate GenBank file 

# Usage
`python3 extract_ti.py [-h] --input INPUT --contiglist CONTIGLIST`

Alternatively, you can submit an SGE job using the wrapper provided in the `SGE_wrapper/run_extract_ti.sh` directory. Use the following command:

`./run_extract_ti.sh path/to/BEAV/output/dir`

# Arguments
- `--input`: A GenBank file generated by BEAV
- `--contiglist`: *.oncogenic_plasmid_final.out.contiglist file that contains all predicted contigs to be oncogenes

# Output
- Ti genes extracted and written in: `*.ti_genes.ffn/fna files`
- Ti gene extraction summary written `in: *.ti_genes.gff3`
- Ti plasmid extracted and saved in: `*.ti_plasmid.gbk`

# To do
- Update bash wrapper script.
- Set up dedicated output directory.
