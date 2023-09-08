#!/usr/bin/python3

import argparse
import beav_oncogenes
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



def get_base_file_name(file_name):
    """
    INPUT: Input GenBank file name or path to the file
    OUTPUT: Base name 
    """
    if '.' in file_name:
        file_name = file_name.split('.')[0]
    
    if '/' in file_name:
        file_name = file_name.split('/')[-1]
    
    return file_name


def get_feature_seq(feature, gbk_record, seq_class, CDS=True):
    '''
    Input: 
        feature: A GenBank sequence feature of interest
        gbk_record: A GenBank record
        seq_class: Type of sequence feature extracting (i.e. T-DNA etc.)
        CDS: If it is a CDS feature or not.
    Output: A description, NA and AA sequences extracted from the record for that feature
    '''
    contig = gbk_record.id
    na_seq = feature.location.extract(gbk_record).seq

    # Get feature locatin
    if feature.strand == -1:
        start = int(feature.location.parts[-1].start)
        end = int(feature.location.parts[0].end)
        strand = '-'
    else:
        start = int(feature.location.parts[0].start)
        end = int(feature.location.parts[-1].end)
        strand = '+'
    
    # Get amino acid sequence and sequence ID
    if CDS != False:
        aa_seq = na_seq.translate(stop_symbol='')
        seq_id = feature.qualifiers['locus_tag'][0]

        if 'gene' in feature.qualifiers.keys():
            gene_id = feature.qualifiers['gene'][0]
        else:
            gene_id = feature.qualifiers['product'][0]
    else:
        seq_id = 'None'
        aa_seq = 'NA'
        gene_id = 'NA'
    

    return (seq_class, contig, seq_id, gene_id, start, end, strand, na_seq, aa_seq)



def extract_ti_features(gbk_file):
    '''
    Input: GenBank file
    Output: List containing corresponding feature sequences and other necessary info
    '''
    feature_list = []
    for gbk_record in SeqIO.parse(gbk_file, "genbank"):
        for feature in gbk_record.features:
            if feature.type == "CDS":
                if 'gene' in feature.qualifiers.keys():
                    if feature.qualifiers['gene'][0] in beav_oncogenes.vir_dict.keys():
                        feature_list.append(get_feature_seq(feature, gbk_record, 'T-DNA Transfer')) #T-DNA transfer
                    elif feature.qualifiers['gene'][0] in beav_oncogenes.oncogene_dict.keys():
                        feature_list.append(get_feature_seq(feature, gbk_record, 'T-DNA')) #T-DNA/Oncogene            
                    elif feature.qualifiers['gene'][0] in beav_oncogenes.rep_dict.keys():
                        feature_list.append(get_feature_seq(feature, gbk_record, 'repABC')) #repABC
                    elif feature.qualifiers['gene'][0] in beav_oncogenes.opine_synth_dict.keys():
                        feature_list.append(get_feature_seq(feature, gbk_record, 'Opine synthase')) #Opine synthase genes
                    elif feature.qualifiers['gene'][0] in beav_oncogenes.opine_tracat_dict.keys():
                        feature_list.append(get_feature_seq(feature, gbk_record, 'Opine transport/catabolism')) #Opine transport/catabolism genes
                    elif feature.qualifiers['gene'][0] in beav_oncogenes.agrocinopine_dict.keys():
                        feature_list.append(get_feature_seq(feature, gbk_record, 'Agrocinopine transport/catabolism')) #Agrocinopine transport/catabolism genes

                else:
                    if feature.qualifiers['product'][0] in beav_oncogenes.oncogene_list:
                        feature_list.append(get_feature_seq(feature, gbk_record, 'T-DNA')) #T-DNA/Oncogene
                    elif feature.qualifiers['product'][0] in beav_oncogenes.rep_dict.values():
                        feature_list.append(get_feature_seq(feature, gbk_record, 'repABC')) #repABC
                    elif feature.qualifiers['product'][0] in beav_oncogenes.opine_synth_dict.values():
                        feature_list.append(get_feature_seq(feature, gbk_record, 'Opine synthase')) #Opine synthase genes
                    elif feature.qualifiers['product'][0] in beav_oncogenes.opine_tracat_dict.values():
                        feature_list.append(get_feature_seq(feature, gbk_record, 'Opine transport/catabolism')) #Opine transport/catabolism genes
                    elif feature.qualifiers['product'][0] in beav_oncogenes.agrocinopine_dict.values():
                        feature_list.append(get_feature_seq(feature, gbk_record, 'Agrocinopine transport/catabolism')) #Agrocinopine transport/catabolism genes
            else:
                if 'note' in feature.qualifiers.keys():
                    if 'origin_of_replication' in ' '.join(feature.qualifiers['note']):
                        feature_list.append(get_feature_seq(feature, gbk_record, 'Origin of replication', CDS=False)) # Origin of replication
                    elif 'virbox' in ' '.join(feature.qualifiers['note']):
                        feature_list.append(get_feature_seq(feature, gbk_record, 'virbox', CDS=False)) # virbox
                    elif 'trabox' in ' '.join(feature.qualifiers['note']):
                        feature_list.append(get_feature_seq(feature, gbk_record, 'trabox', CDS=False)) # trabox
                    elif 'T-DNA_right_border' in ' '.join(feature.qualifiers['note']):
                        feature_list.append(get_feature_seq(feature, gbk_record, 'T-DNA right border', CDS=False)) # t-dna right border
                    elif 'T-DNA_left_border' in feature.qualifiers['note']:
                        feature_list.append(get_feature_seq(feature, gbk_record, 'T-DNA left border', CDS=False)) # t-dna left border

    return feature_list


def make_SeqRecord(feature_entry, seq_type='NA'):
    '''
    Input: Single list containing feature and sequences
    Output: A SeqRecord object
    '''
    if seq_type == 'NA':
        seq = feature_entry[7] # NA seq record
    elif seq_type == 'AA':
        seq = feature_entry[8] # AA seq record

    seq_record = SeqRecord(
            seq, # NA Seq Object
            id = '', # Locus Tag
            description = f"id={feature_entry[2]} | contig={feature_entry[1]} | type={feature_entry[0]} | gene={feature_entry[3]} | start={feature_entry[4]} | end={feature_entry[5]} | strand={feature_entry[6]}"
        )

    return seq_record    


def make_ti_genes_fasta(feature_list, input):
    '''
    Input: List containing features and sequences from different contigs/GenBank records
    Output: Write fasta file containing NA and AA sequences
    '''
    na_seqs = []
    aa_seqs = []

    for f in feature_list:
        if f[2] != 'None':
            na_seqs.append(make_SeqRecord(f, 'NA'))
            aa_seqs.append(make_SeqRecord(f, 'AA'))

    SeqIO.write(na_seqs, f"{get_base_file_name(input)}.ti_genes.ffn", 'fasta')
    SeqIO.write(aa_seqs, f"{get_base_file_name(input)}.ti_genes.faa", 'fasta')


    return f"Ti genes extracted and written in: f{get_base_file_name(input)}.ti_genes.ffn/fna files" 


def make_gff3(features_list, input):
    '''
    Input: Feature list from get_feature_seq function
    Output: GFF3 summarizing which contig has which genes  
    '''

    with open(f"{get_base_file_name(input)}.ti_genes.gff3", 'w') as gff:
        gff.writelines('##gff-version 3\n')
        for f in features_list:
            if f[2] != 'None':
                f_type = 'CDS' 
            else:
                f_type = 'misc' 
            # gff columns: seqid, source, type, start, end, score, strand, phase, attributes
            line = f"{f[1]}\tBEAV\t{f_type}\t{f[4]}\t{f[5]}\t.\t{f[6]}\t.\ttype={f[0]};locus_tag={f[2]};gene={f[3]}"
            gff.writelines(f"{line}\n")

    return f'Ti gene extraction summary written in: {get_base_file_name(input)}.ti_genes.gff3'


def extract_ti_plasmid(gbk_file, contiglist):
    '''
    Output: GenBank file containing ti plasmid contig(s)
    '''

    # Find the contigs that contain Ti Plasmid
    with open(contiglist, 'r') as f:
        raw = f.readlines()
    
    contigs = [i.split()[0] for i in raw]

    # Filter beav gbk file to have only oncogenic contigs

    records = []

    for gbk_record in SeqIO.parse(gbk_file, "genbank"):
        if gbk_record.id in contigs:
            records.append(gbk_record)
    
    # Write oncogenic plasmid record
    SeqIO.write(records, f"{get_base_file_name(gbk_file)}.ti_plasmid.gbk", 'genbank')

    return  f"Ti plasmid extracted and saved in: {get_base_file_name(gbk_file)}.ti_plasmid.gbk"


def main():
    # Define CLI parser argument
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', '-i', type=str, required=True)
    parser.add_argument('--contiglist', type=str, required=True)

    # Parse the arguments
    args = parser.parse_args()

    # Extract Ti associated genes and write to fna/faa/gff3 files
    ti_features = extract_ti_features(args.input)
    print(make_ti_genes_fasta(ti_features, args.input))
    print(make_gff3(ti_features, args.input))

    # Extract whole ti contigs
    print(extract_ti_plasmid(args.input, args.contiglist))


if __name__ == "__main__":
    main()
