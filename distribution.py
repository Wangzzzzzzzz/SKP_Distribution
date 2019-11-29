from blastp import obtain_seq, write_to_fasta
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from io import StringIO
import os, math
import numpy as np
import pandas as pd 

def run_blast(query, subject, evalue):
    """
    run blast with query against subject, keeping result larger than evalue
    """
    output = NcbiblastpCommandline(
        query=query, subject=subject, outfmt=5, evalue=evalue)()[0]
    # read all parser result
    blast_result_record = list(NCBIXML.parse(StringIO(output)))
    return blast_result_record

def read_fasta(fasta_file):
    return [record for record in SeqIO.parse(fasta_file,"fasta")]

def organize_blast_res(blast_res, query_file, subject_file, **kwargs):
    """
    organize the blast result into a table with the 
    following format:
    subject_seq, query_seq, score, bit_score, E-value, log-E-value, per_identity

    Pass the max/min value of the score, bit_score, E-value, log-E-value, per_identity
    into the function if neccessary.

    If not given the following default will be used:
    min_score = 0
    min_bit_score = 0
    max_E_value = Inf
    max_log_E_value = Inf
    min_per_identity = 0
    """
    query_file = './fastaDB/' + query_file
    subject_file = './fastaDB/' + subject_file
    user_def_max_min = set(kwargs.keys())
    blast_rs_dict = {}
    for query_seq in blast_res:
        query_seq_name = query_seq.query[0:6]
        blast_rs_dict[query_seq_name] = dict()
        # go through all alignments of all one sequence in query
        for alignment in query_seq.alignments:
            # subject seq (aligned sequence in database)'s name
            subjt_seq_name = alignment.title[0:6]

            # that contains significant alignment, usually there is only one)
            hsp_list = []
            for hsp in alignment.hsps:
                hsp_discriptor = {
                    "score": hsp.score,
                    "E_Value": hsp.expect,
                    "Log_E_Value": float('-inf') if hsp.expect == 0 else np.log(hsp.expect),
                    "Bit_Score": hsp.bits,
                    "Perc_Identity": hsp.identities/hsp.align_length*100,
                    "Num_of_Identity": hsp.identities,
                    "Subject_Name": subjt_seq_name,
                    "Query_Name": query_seq_name
                }
                hsp_list.append(hsp_discriptor)
            # set the result to be the best match based on bit score
            blast_rs_dict[query_seq_name][subjt_seq_name] = max(hsp_list, key=lambda x: x['Bit_Score'])
    query_name = sorted([seq.name for seq in read_fasta(query_file)])
    subject_name = sorted([seq.name for seq in read_fasta(subject_file)])
    query_check_set = set(blast_rs_dict.keys())
    organized_res = []
    for query in query_name:
        if query in query_check_set:
            subj_check_set = set(blast_rs_dict[query].keys())
            for subjt in subject_name:
                if subjt in subj_check_set:
                    discriptor = blast_rs_dict[query][subjt]
                    organized_res.append((
                        subjt,
                        query,
                        discriptor['score'],
                        discriptor['Bit_Score'],
                        discriptor['E_Value'],
                        discriptor['Log_E_Value'],
                        discriptor['Perc_Identity']
                    ))
                else:
                    organized_res.append((
                        subjt,
                        query,
                        kwargs["min_score"] if "min_score" in user_def_max_min else 0,
                        kwargs["min_bit_score"] if "min_bit_score" in user_def_max_min else 0,
                        kwargs["max_E_value"] if "max_E_value" in user_def_max_min else float('Inf'),
                        kwargs["max_log_E_value"] if "max_log_E_value" in user_def_max_min else float('Inf'),
                        kwargs["min_per_identity"] if "min_per_identity" in user_def_max_min else 0
                    ))
        else:
            for subjt in subject_name:
                organized_res.append((
                    subjt,
                    query,
                    kwargs["min_score"] if "min_score" in user_def_max_min else 0,
                    kwargs["min_bit_score"] if "min_bit_score" in user_def_max_min else 0,
                    kwargs["max_E_value"] if "max_E_value" in user_def_max_min else float('Inf'),
                    kwargs["max_log_E_value"] if "max_log_E_value" in user_def_max_min else float('Inf'),
                    kwargs["min_per_identity"] if "min_per_identity" in user_def_max_min else 0
                ))
    return organized_res
        



def main():
    if not os.path.exists('./dataFile'):
        os.mkdir('./dataFile')
    if not os.path.exists('./fastaDB'):
        os.mkdir('./fastaDB')
    if not os.path.exists('./output'):
        os.mkdir('./output')

    v1_v2_col_names = ['SKP_v1_wild', 'SKP_v2_wild', 'score',
                       'bit_score', 'E_value', 'log_E_value', 'per_identity']
    identity_col_names = ['wild_1', 'wild_2', 'score',
                          'bit_score', 'E_value', 'log_E_value', 'per_identity']
    
    # process the file for skempi v1 and skempi v2
    skempi_v1 = obtain_seq('SKP1102s.ddg.txt', 'SKP1102s.seq.txt')
    write_to_fasta(skempi_v1, './fastaDB/skempi_v1.fasta')
    skempi_v2 = obtain_seq('3G_S487_test_dataset_top1.txt','skempi_v2.singlemut.mut4.seq.txt')
    write_to_fasta(skempi_v2, './fastaDB/skempi_v2.fasta')
    
    # run blast
    skempi_v1_v1 = run_blast(query='./fastaDB/skempi_v1.fasta',
                             subject='./fastaDB/skempi_v1.fasta',
                             evalue=1e20)
    skempi_v2_v2 = run_blast(query='./fastaDB/skempi_v2.fasta',
                             subject='./fastaDB/skempi_v2.fasta',
                             evalue=1e20)
    skempi_v2_v1 = run_blast(query='./fastaDB/skempi_v2.fasta',
                             subject='./fastaDB/skempi_v1.fasta',
                             evalue=1e20)

    # organize results
    org_v2_v1 = organize_blast_res(blast_res=skempi_v2_v1,
                                   query_file='skempi_v2.fasta',
                                   subject_file='skempi_v1.fasta') 
    org_v1_v1 = organize_blast_res(blast_res=skempi_v1_v1,
                                   query_file='skempi_v1.fasta',
                                   subject_file='skempi_v1.fasta')
    org_v2_v2 = organize_blast_res(blast_res=skempi_v2_v2,
                                   query_file='skempi_v2.fasta',
                                   subject_file='skempi_v2.fasta')
    
    # output to datasets
    dataset_v2_v1 = pd.DataFrame(org_v2_v1,columns=v1_v2_col_names)
    dataset_v1_v1 = pd.DataFrame(org_v1_v1,columns=identity_col_names)
    dataset_v2_v2 = pd.DataFrame(org_v2_v2,columns=identity_col_names)
    dataset_v2_v1.to_csv('output/SKP_V1_V2.tsv','\t', index=False)
    dataset_v1_v1.to_csv('output/SKP_V1_V1.tsv', '\t', index=False)
    dataset_v2_v2.to_csv('output/SKP_V2_V2.tsv', '\t', index=False)

if __name__ == '__main__':
    main()
