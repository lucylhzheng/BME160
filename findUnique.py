#!/usr/bin/env python3
#Name: Lucy Zheng (lzheng20)

'''
The findUnique module takes in fasta file of tRNA sequences and generates the essential sets for each tRNA.

Expected input at the command line:

    python3 findUnique.py <bos-tRNA-7.fa> output.txt

Expected output:

tRNA | Ala | UGC | Bos taurus | mitochondrial
GAGGAUUU"LCUUAAUUAAAGULGPUGAUUUGCAUPCAAUUGAUGUAAGGUGPAGUCUUGCAAUCCUUACCA
GAGGA
..GGAU
...GAUUU"
......UU"LC
.........LCUUAAU
...........UUAAUUA
..............AUUAAA
...............UUAAAG
...................AGUL
.....................ULG
......................LGP
.......................GPUGA
........................PUGAU
...........................AUUUG
............................UUUGC
..............................UGCAU
.................................AUPC
...................................PCAA
....................................CAAUUG
.....................................AAUUGA
......................................AUUGAUG
.........................................GAUGUA
...........................................UGUAAGG
...............................................AGGUG
...................................................GPA
....................................................PAGU
.....................................................AGUCUU
......................................................GUCUUG
........................................................CUUGC
.........................................................UUGCAA
...........................................................GCAAU
.............................................................AAUC
..............................................................AUCCUUA
................................................................CCUUAC
(continues with next tRNA header)

'''
class trnaClass :
    '''
    Generates the header, sequence, and a list of essential and unique elements for each tRNA.
    '''
    tRNAs = []
 #stores each tRNA header and its sequence 
    def __init__(self,header,tseq):
        '''
        Stores lists used for each class method and filter sequences.

        Args:
            header (str): tRNA headers are stored in the trnaClass.tRNAs list.
            tseq (str): tRNA sequences are stored in the trnaClass.tRNAs list.
        '''
        seq = tseq.replace("-","").replace(".","").replace("_","")
 # takes out the [-._] characters
        trnaClass.tRNAs.append((header,seq))
 # appending header and seq to trnaClass.tRNAs list
        self.pset_list = []
        self.unique = []
        self.essentials = []
        self.index = []
        
    def Addsequence (self):
        '''
        Generates powersets for each tRNA sequence

        Returns:
            A list of powersets for each tRNA in self.pset_list.
        '''
        for header, seq in trnaClass.tRNAs:
            power_set = set()
            for start in range(0,len(seq)):
                for end in range(start,len(seq)):
                    substrs = seq[start:end+1]
                    power_set.add(substrs)
 #adding each set to power_set
            self.pset_list.append(power_set)   
    def Uniquebuilder(self):
        '''
        Identifies and creates sets of unique elements for each tRNA sequence.

        Returns:
            A list with sets of unique elements for each tRNA.
        '''
        for index, main_tRNA in enumerate(self.pset_list):
            other_tRNAs = set()
            for index_, subset in enumerate(self.pset_list):
                if index != index_:
                    other_tRNAs = other_tRNAs.union(subset)
 # Appending all the set in the powerset that is not within the current tRNA
            main_unique = main_tRNA.difference(other_tRNAs)
            self.unique.append(main_unique)
    def Essentialsbuilder(self):
        '''
        Sort through the unique sets of each tRNA to find the essential elements.

        Returns:
            Unique sets of each tRNA in the self.essentials list.
        '''
        for main_set in self.unique:
            small_set = set()
            for seq in main_set:
                if seq[1:] in main_set:
 
                    small_set.add(seq)
                elif seq[:-1] in main_set:
                    small_set.add(seq)
            small_set = main_set - small_set
 #uses the equation "essentials = unique - nonessentials"
            self.essentials.append(small_set)
    def PrintEssentials(self):
        '''
        Finds index for each unique element and outputs in the desired format with the tRNA headers and sequences.

        Returns:
            Each tRNA has their own output of its header, sequence, and essential elements.
        
        '''
        essentials = []
        for index, ess in enumerate(self.essentials):
            var = []
            header, seq = trnaClass.tRNAs[index]
            for substrings in ess:
                test = 0
                position = 0
                while(test != -1):
 # looking for the index of each essential element
                    test = seq[position:].find(substrings)
                    if test != -1:
 # testing for an element that shows up more than one time in the sequence
                        var.append((substrings, position+test))
                        position = position + test + 1
            essentials.append(sorted(var, key = lambda element:element[1])) # sorts the essential elements based on its index number
            output = zip(trnaClass.tRNAs,essentials)
            final = list(output)
        for seq,ele in sorted(final):
            print(seq[0])
            print(seq[1])
            for x,y in ele:
                dots = "".join(['.'] * y)
                print("{0}{1}".format(dots,x))
import sys
class FastAreader :
    '''
    Parses through the header and sequence of each tRNA in the fasta file.
    '''
    
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        '''opens the file inserted from commandline '''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' identifies and extracts the header and its sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                if not line: # we are at EOF
                    return header, sequence
                line = fileH.readline()
            header = line[2:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header,sequence
                    header = line[2:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

########################################################################
# Main
# Here is the main program
########################################################################
def main(inCL=None):
    '''
    Take each tRNA sequence and generate its essential sets.

    Return:
        All 22 tRNA has their output of each tRNA header, sequence, and essential sets.
    '''
    myReader = FastAreader()
    for header, seq in myReader.readFasta():
 # passing each header and sequence to trnaClass
        mytrnaClass = trnaClass(header,seq)
    mytrnaClass.Addsequence()
    mytrnaClass.Uniquebuilder()
    mytrnaClass.Essentialsbuilder()
    mytrnaClass.PrintEssentials()
    
if __name__ == "__main__":
    main() 
