#!/usr/bin/python
import argparse
import sys
from Bio import SeqIO
import numpy as np
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Construction of a distance matrix between each pair of sequences in FASTA file.')
    parser.add_argument('-i', '--input', help='Input FASTA file', metavar='FILE', required=True)
    parser.add_argument('-t', '--threads', help='Number of threads', metavar='Int', type=int, default=1)
    parser.add_argument('-o', '--output', help='The resulting filename', metavar='FILE', default='/home/izmaylova/Documents/Matrix.odt')
    args = parser.parse_args()
    # Creation of dictionary for penalties computation for nucleotide substitutions and gaps.
    subs_matrix = {}
    subs_matrix['A', 'A']=10
    subs_matrix['A', 'G']=-1
    subs_matrix['A', 'C']=-3
    subs_matrix['A', 'T']=-4
    subs_matrix['G', 'A']=-1
    subs_matrix['G', 'G']=7
    subs_matrix['G', 'C']=-5
    subs_matrix['G', 'T']=-3
    subs_matrix['C', 'A']=-3
    subs_matrix['C', 'G']=-5
    subs_matrix['C', 'C']=9
    subs_matrix['C', 'T']=0
    subs_matrix['T', 'A']=-4
    subs_matrix['T', 'G']=-3
    subs_matrix['T', 'C']=0
    subs_matrix['T', 'T']=8
    all_letters = set( [key[0] for key in subs_matrix.keys()] )
    subs_matrix.update(  {(letter,'-'):-5 for letter in all_letters} )
    subs_matrix.update(  {('-',letter):-5 for letter in all_letters} )
    # Function for score calculation between two compared sequences.
    def score(seq1, seq2):
        n, m = len(seq1), len(seq2)
        Wmax = np.zeros( (n+1, m+1) )
        for i in range(n+1):
            for j in range(m+1):
                if 0 <= i-1 < n and 0 <= j-1 < m:
                    a1 = Wmax[i, j-1]+subs_matrix['-', seq2[j-1]]
                    a2 = Wmax[i-1, j]+subs_matrix[seq1[i-1], '-']
                    a3 = Wmax[i-1, j-1]+subs_matrix[seq1[i-1], seq2[j-1]]
                    Wmax[i,j] = max([a1, a2, a3])
                elif i-1 == -1 and 0 <= j-1 < m:
                    Wmax[i,j] =  Wmax[i, j-1]+subs_matrix['-', seq2[j-1]]
                elif j-1 == -1 and 0 <= i-1 < n:
                    Wmax[i,j] = Wmax[i-1, j]+subs_matrix[seq1[i-1], '-']
        return Wmax[n, m]
    # Recording of FASTA file sequences to the list.
    records = list(SeqIO.parse(args.input, "fasta"))
    dim = len(records)
    # Creation of a score matrix for distance computation between each pair of sequences in FASTA file.
    Wres = np.zeros( (dim, dim) )
    for i in range(dim):
        for j in range(dim):
            Wres[i,j] = score(str(records[i].seq).upper(), str(records[j].seq).upper())
    sys.stdout = open(args.output, 'w')
    print(Wres)
    sys.stdout.close()

