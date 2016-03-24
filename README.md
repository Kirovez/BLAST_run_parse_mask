# BLAST_run_parse_mask
Python 3.4 script to run BLASTn and create  1) fasta file contained query sequences showed similarity to sbjct; 2) BLAST_MASK file with query sequences with masked 'n' regions participated in alignment;  3) fasta file of queries which did not show similarity; 4) HIT fasta file contained parts of sequences participated in alignment. 
__author__ = 'ikirov'
from Bio import SeqIO
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SearchIO

list=[]
in_file =  ''
DB_file = ''
Similarity = 90
Evalue = 1e-5
no_hits_ind = 0
query_list = []
Coverage = 90
query_indexed_db  = SeqIO.index(in_file,'fasta')
#sbjct_indexed_db  = SeqIO.index(DB_file,'fasta')
print('Number of query sequences:' + str(len(query_indexed_db)))

file_name = '{0}_vs_{1}'.format(in_file,DB_file)
with open(file_name, 'w') as file:
    print('File opening')
cmd = r'C:\Users\ikirov.CLO\Desktop\blast-2.2.31+\bin\makeblastdb.exe -in %s -dbtype nucl' % DB_file
os.system(cmd)
proga = r'C:\Users\ikirov.CLO\Desktop\blast-2.2.31+\bin\blastn.exe'
blast = NcbiblastnCommandline(proga, query=in_file, db=DB_file, out=file_name, outfmt=5, word_size = 12, evalue = 0.001, strand = 'plus')#strand = 'plus'
stdout, stderr = blast()

list_of_hit_numbers = []
list_sbjct = []

i=0
g=0
cal = 0
index=0
query_count = 0
number_of_quries_in_fasta_file =0
list_of_printed=[]
sbjct_list = []

with open('HIT', 'w') as hit, open ('BLAST_MASK', 'w') as BLAST_MASK, open('query_in_fasta','w') as query_in_fasta, open('NO_hits', 'w') as no_hits:
    file2 = open(file_name)
    s1 = SearchIO.parse(file2, 'blast-xml')
    count_hits = 0
    for recor in s1:
        number_of_hits=0
        query_count+=1
        print(query_count)
        for HSP in recor:
            #sbj = str(sbjct_indexed_db[HSP.id].seq)
            sequence = str(query_indexed_db[recor.id].seq)
            indicator = 0
            hs = HSP.hsps

            for u in hs:
                if u.evalue<=Evalue and u.ident_num*100/len(u.query) >= Similarity:
                    sequence = sequence.replace(sequence[u.query_range[0]:u.query_range[1]], 'n' * len(u.query))
                    indicator += 1
                    g+=1
                    count_hits+=1
                    #if len(u.query)>=(Coverage*len(sequence))/100 or len(u.hit)>=(Coverage*len(sbj))/100:
                    SeqIO.write(u.query,hit,'fasta')
                    SeqIO.write(u.hit,hit,'fasta')
                    cal+=1
            coverage = sequence.count('n')*100/len(sequence)
            #Coverage calculation based on the masking results
            if coverage>=Coverage and indicator >= 1:
                number_of_hits+=1
                BLAST_MASK.write(('>' +recor.id+'_'+str(HSP.id) + '\n' + str(sequence) + '\n'))
                if HSP.id != recor.id and HSP.id not in list_sbjct:
                        list_sbjct.append(HSP.id)
                if recor.id not in list_of_printed:
                    number_of_quries_in_fasta_file +=1
                    SeqIO.write(query_indexed_db[recor.id],query_in_fasta,'fasta')
                    list_of_printed.append(recor.id)
                    query_list.append(recor.id)

        if number_of_hits==0:
            no_hits_ind+=1

    # NO hits file writing

    for seq in SeqIO.parse(in_file,'fasta'):
        if seq.id not in query_list:
            SeqIO.write(seq,no_hits,'fasta')


print('Hit pairs in HIT file: ' + str(cal))
print('Number of query sequences in fasta file: '+str(number_of_quries_in_fasta_file))
print('Number of queries with no hits: ' + str(no_hits_ind))
print('Number of different sbjct sequences: ' + str(len(list_sbjct)))
print(query_list)


