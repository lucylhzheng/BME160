#!/usr/bin/env python3 
# Name: Lucy Zheng (lzheng20)

'''
Sequence information is inserted and gets converted to the respective dictionaries of RNA,DNA, and amino acids which outputs the 
the correct conversion based on the data given

examples:
Input: “ATG” 
output: ATG = MET

Input: “UAG” 
output: UAG = ---

Input: “E” 
output: E = GLU

Iniput: “Asp”
output: ASP = D

'''
short_AA = {
            'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
            } 

long_AA = {value:key for key,value in short_AA.items()} #reverse key and value mappings for one character code

rnaCodonTable = {
# Second Base
# U             C             A             G
#U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
#C 
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
#A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
'AUG': 'MET', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
#G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}
dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()} #replaces "U" for "T" after mapping the key with the value

def main():
    ''' sequence information is inputted and go through IF statements to undergo the targeted converted '''
    code = input("Enter an amino acid:")
    if len(code) == 1: #one character code goes through the amino acid converted
        print(code + " = " + (long_AA.get(code.upper(),"unknown")))
    elif "T" in code.upper(): #string with "T" will require mapping from the DNA converter to use the RNA dictionary
        print(code + " = " + dnaCodonTable.get(code.upper(),"unknown"))
    elif "U" in code: #string with "U" needs mapping from the RNA dictionary
        print(code + " = " + rnaCodonTable.get(code.upper(),"unknown"))
    else:
        print(code + " = " + short_AA.get(code.upper(),"unknown"))

    
    
main()
