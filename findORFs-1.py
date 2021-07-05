#!/usr/bin/env python3
# Name: Lucy Zheng (lzheng20)
# Group Members: Christopher Tran (Chlotran)


########################################################################
# Design:
# The three classes (orfFinder, FastAreader, and NucParams are imported from sequenceAnalysis module at the beginning of the program.
# Design of the findORF class includes three class methods (reverse_comp, top_strand, and bottom_strand).
# Reverse_comp class method is used to find the complentary strand of the given DNA sequence using the replace command and reverse indexing.
# The top_strand and bottom_strand find the targetted gene within each frame by using a for loop that stores all the start codon found in a list.
# To target the boundary conditions, if statements are used to check for two conditions where there are no start codons and with frame that does not have either start and stop codon.
########################################################################
from sequenceAnalysis import OrfFinder, FastAreader
'''
The findORFs module takes in fasta files and specific parameters using CommandLine class and generates different open reading frames in the imported OrfFinder class and FastAreader.

Expected input at the command line:
python3 findORFs <lab5test.fa> output.txt -lG=False -mG=9

Expected output:

test

+1     1..    9     9
+3     1..    9     9
-1     1..    9     9
-2     1..    9     9
-3     1..    9     9
test2
+1     1..   10    10
-1     1..   10    10
-2     1..   10    10
-3     1..   10    10
+2     2..   10     9
test3
+2     1..   11    11
-1     1..   11    11
-2     1..   11    11
-3     1..   11    11
+3     3..   11     9
test-1
+1     1..    9     9
+2     1..    9     9
+3     1..    9     9
-1     1..    9     9
-3     1..    9     9
test-2
+1     1..   10    10
+2     1..   10    10
+3     1..   10    10
-1     1..   10    10
-2     1..    9     9
test-3
+1     1..   11    11
+2     1..   11    11
+3     1..   11    11
-2     1..   11    11
-3     1..    9     9
test1A
+3     1..   10    10
-1     1..   10    10
-2     1..   10    10
-3     1..   10    10
+1     1..    9     9
test2A
+1     1..   11    11
-1     1..   11    11
-2     1..   11    11
-3     1..   11    11
+2     2..   10     9
test3A
+2     1..   12    12
-1     1..   12    12
-2     1..   12    12
-3     1..   12    12
+3     3..   11     9
test-1A
+1     1..   10    10
+2     1..   10    10
+3     1..   10    10
-3     1..   10    10
-1     2..   10     9
test-2A
+1     1..   12    12
+2     1..   12    12
+3     1..   12    12
-1     1..   12    12
-2     3..   11     9
test-3A
+1     1..   13    13
+2     1..   13    13
+3     1..   13    13
-2     1..   13    13
-3     3..   11     9
    
'''
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, # default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        #self.parser.add_argument('inFile', action = 'store', help='input file name')
        #self.parser.add_argument('outFile', action = 'store', help='output file name') 
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', help='start Codon') # allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


def main(inCL=None):
    '''
    Input given DNA sequences into orfFinder class and formats the output given.

    Args:
        inCL: checks for input given from the command line.
        
    Return:
        No output if inCL is set to none and open reading frames are outputted is inCL is not set to none.
    '''
    if inCL is None:
        myCommandLine = CommandLine()
                        
    else :
        myCommandLine = CommandLine(inCL)
    
    longest_gene = True
    if myCommandLine.args.longestGene == "False":
        longest_gene = False # setting condition for the Longest gene
    myReader = FastAreader()
    myOrfFinder = OrfFinder(lg=longest_gene,mg=myCommandLine.args.minGene,start=myCommandLine.args.start,stop=myCommandLine.args.stop)
    for head, seq in myReader.readFasta():
        print(head)
        topSeqs = myOrfFinder.topStrand(seq)
        reverseSeq = myOrfFinder.reverseComp(seq)
        bottomSeqs = myOrfFinder.bottomStrand(reverseSeq)
        all_orfs = []
        all_orfs.extend(topSeqs)
        all_orfs.extend(bottomSeqs)
        def sortOrfs(nucs):
            '''
            Sorts all the output based on the largest to smallest of the gene length and starting position.
        
            Args:
                orfs: all the orfs outputted from the OrfFinder class
            '''
            lens= nucs[2] - nucs[1] # descending order of the ORFs
            return (-1*lens,nucs[1])
        all_orfs.sort(key=sortOrfs) 
        for x in all_orfs:
            fr = x[0]
            sta = x[1]
            stp = x[2]
            l = x[3]
            if l < myCommandLine.args.minGene: #accounting for minimum gene value and extra credit of varying min gene sizes
                continue
            if fr < 0:
                print("{:d} {:>5d}..{:>5d} {:>5d}".format(fr,sta,stp,l))
            else:
                print("{:+d} {:>5d}..{:>5d} {:>5d}".format(fr,sta,stp,l))

    #sort by size then starting position
    
if __name__ == "__main__":
    main()  # delete the list when you want to run with STDIN
    
