#!/bin/bash
infile=$1

eval "$(/nfs7/BPP/Weisberg_Lab/lab_members/rahman/opt/miniconda3/condabin/conda shell.bash hook)"

conda activate py38

python3 extract_ti.py --input ${infile}/*gbk --contiglist ${infile}/*oncogenic_plasmid_final.out.contiglist

conda deactivate
