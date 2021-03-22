# !/usr/bin/env python3
# Name: Lucy Zheng (lzheng20)
# Group Members: Christopher Tran (chlotran)
'''
This program takes the given data that is altered from FastAreader and calls the
specific classes and methods within SequenceAnalysis module for proper formatting
of the desired output.
'''


from sequenceAnalysis import NucParams, FastAreader, ProteinParam # imports specific class methods in module

def main ():
    '''
    This function takes the given DNA sequences and format the ouputs from class methods in sequenceAnalysis module
    which is called and sorts all the information in the output format desired.
    '''
    
    myReader = FastAreader() # make sure to change this to use stdin
    myNuc = NucParams()
    seq = myReader.readFasta()
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)

    total_nucs = sum(myNuc.nucComp.values())
    seq_len = total_nucs/1000000
    print("sequence length = {:.2f} Mb".format(seq_len)) # outputs the total nucleotides sequence in Mb
    print() # empty line space
     
    
    g_nucs = myNuc.nucComp["G"]
    c_nucs = myNuc.nucComp["C"]
    gc_nucs = g_nucs + c_nucs
    gc_comp = (gc_nucs/total_nucs) * 100
    print("GC content = {:.1f}%".format(gc_comp)) # gives the GC content in the given sequences as a percentage
    print()
    
    
    # sort codons in alpha order, by Amino Acid
    # calculate relative codon usage for each codon and print
    for codon, aa in sorted(myNuc.rnaCodonTable.items(),key=lambda x: x[1]): # sorts the amino acids in alphabetical order
        aa_counts = myNuc.aaComp[aa]
        codon_counts = myNuc.codonComp[codon]
        val = (codon_counts/aa_counts)

        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, aa, val*100, myNuc.codonComp[codon]))

if __name__ == "__main__":
    main()

