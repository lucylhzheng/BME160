# !/usr/bin/env python3
# Name: Lucy Zheng (lzheng20)
# Group Members: Christopher Tran (chlotran)
'''
SequenceAnalysis is a module that contains three classes which are NucParams, FastAreader, and ProteinParams
which can be called from the program titled GenomeAnalyzer. All the calculations of amino acids, nucleotides,
and codons from the data inputted will be done here only if it is called from GenomeAnalyzer and will produce
the output.
'''

class NucParams:
    '''
    NucParams is where the calculations and data storing occurs which are
    categorized into amino acids, codons, and nucleotides.
    '''
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G',  # GxG
    }
    
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
    
    def __init__ (self, inString=''):
      '''
      This is where all the three essential dictionaries of nucleotides,
      amino acids, and codons are created, stored, and used.
      '''
      
      self.inString = inString.upper() # capitalized all characters in the input string
    
      aaComp = {n:0 for n in "ACDERGHILKMNPQRSTVYW-"}
      self.aaComp = aaComp # dictionary for amino acids

      nucComp = {}
      self.nucComp = nucComp # dictionary for nucleotides
    
      codonComp = {
      # U
      'UUU': 0 , 'UCU': 0, 'UAU': 0, 'UGU': 0,  
      'UUC': 0 , 'UCC': 0, 'UAC': 0, 'UGC': 0,  
      'UUA': 0 , 'UCA': 0, 'UAA': 0, 'UGA': 0,  
      'UUG': 0 , 'UCG': 0, 'UAG': 0, 'UGG': 0,  
      # C
      'CUU': 0, 'CCU': 0, 'CAU': 0, 'CGU': 0,  
      'CUC': 0, 'CCC': 0, 'CAC': 0, 'CGC': 0,  
      'CUA': 0, 'CCA': 0, 'CAA': 0, 'CGA': 0, 
      'CUG': 0, 'CCG': 0, 'CAG': 0, 'CGG': 0,  
      # A
      'AUU': 0, 'ACU':  0, 'AAU': 0, 'AGU': 0,  
      'AUC': 0, 'ACC':  0, 'AAC': 0, 'AGC': 0,  
      'AUA': 0, 'ACA':  0, 'AAA': 0, 'AGA': 0,  
      'AUG': 0, 'ACG':  0, 'AAG': 0, 'AGG': 0, 
      # G
      'GUU': 0, 'GCU': 0, 'GAU': 0, 'GGU': 0,  
      'GUC': 0, 'GCC': 0, 'GAC': 0, 'GGC': 0,  
      'GUA': 0, 'GCA': 0, 'GAA': 0, 'GGA': 0,  
      'GUG': 0, 'GCG': 0, 'GAG': 0, 'GGG': 0,
      } 
      self.codonComp = codonComp # dictionary for codon
    
    def addSequence (self, inSeq):
      '''
      This class method updates all three dictionary from
      the __init__ method from the given data through specific
      data conversions.
      '''

      InSeq = inSeq.upper()

      # aaComposition
      rna_string = ['A','C','G','T','U']
      nuc_list = []
      for n in InSeq:
          if n in rna_string:
              nuc_list.append(n)
      nuc_string = ''.join([str(elem) for elem in nuc_list])
      rna = nuc_string.replace("T","U")
      for codon in range(0,len(rna),3):
          x = self.rnaCodonTable.get(rna[codon:codon+3],'')
          self.aaComp[x] = self.aaComp.get(x,0)+1
          

      # nucComposition
      string = ['A','C','G','T','U','N']
      for n in InSeq:
          if n in string:
            self.nucComp[n] = self.nucComp.get(n,0)+1

      # codonComposition
      rna = InSeq.replace("T","U").replace(' ','')
      for codon in range(0,len(rna),3):
          rna_codons = rna[codon:codon+3]
          if "N" in rna_codons and len(rna_codons)==3:
              continue
          elif rna_codons in self.codonComp:
              self.codonComp[rna_codons] += 1
          else:
              self.codonComp[rna_codons] += 0


    def aaComposition(self):
      '''
      This class method ouputs the amino acid dictionary.
      '''
      return self.aaComp
    
    def nucComposition(self):
      '''
      This class method outputs the nucleotides dictionary.
      '''
      return self.nucComp
    def codonComposition(self):
      '''
      This class method outputs the codons dictionary.
      '''
      return self.codonComp
    def nucCount(self):
      '''
      This class method outputs the total sum of the nucleotides
      within the given data.
      '''
      return sum(self.nucComp.values())
		

import sys
class FastAreader:
  '''This class reads the data from the given file
  and extract the dna sequences without the given headers.'''
  
  def __init__ (self, fname=''):
    '''contructor: saves attribute fname '''
    self.fname = fname
  def doOpen (self):
    ''' Handle file opens, allowing STDIN.'''
    if self.fname == '':
      return sys.stdin
    else:
      return open(self.fname)
  def readFasta (self):
    ''' Read an entire FastA record and return the sequence header/sequence'''
    header = ''
    sequence = ''
    with self.doOpen() as fileH:
            
      header = ''
      sequence = ''
            
      # skip to first fasta header
      line = fileH.readline()
      while not line.startswith('>') :
       line = fileH.readline()
      header = line[1:].rstrip()
      for line in fileH:
        if line.startswith ('>'):
            yield header,sequence
            header = line[1:].rstrip()
            sequence = ''
        else :
            sequence += ''.join(line.rstrip().split()).upper()
    yield header,sequence


class ProteinParam:
  '''This class takes a string of amino acids and calculates
  the physical-chemical properties related to the given protein sequence.'''
  aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.20, 'Y': 181.189
        }
  mwH2O = 18.015
  aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
  aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
  aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
  aaNterm = 9.69
  aaCterm = 2.34

  def __init__ (self, protein):
    '''Computes and save the aaComposition Dictionary of the inputted protein sequence.
    Protein is also saved an attribute'''
    self.protein = protein
    aaComp = {n:0 for n in "ACDERGHILKMNPQRSTVYW"}
    for n in protein:
      aaComp[n] = protein.count(n)
      self.aaComp = aaComp
  
  def aaCount (self):
    '''Counts every valid amino acid from the protein sequence and outputs the integer of the total count'''
    cap_proteins = self.protein.upper()
    substrings = ["A", "C", "D", "E", "F", "G", "H", "I", "L", "K", "M", "N", "P", "Q" , "R" , "S", "T", "V", "Y", "W"]
    count = 0 
    for protein in cap_proteins:
      if protein in substrings: # If statement is used after the for loop to count the amino acid that is also listed in the substrings
        count += 1 # adds 1 every time, there is a valid amino acid within the protein sequence
      else:
        count += 0 # If amino acid does not match elements in substrings, it will not be counted
        break
        return count
          
  def pI (self):
    '''finds and outputs the specific pH that is close to the neutral net charge (0)'''
    lowest_ph = 0.0
    lowest_charge = abs(self._charge_(0.0))
    for x in range(0,1401): # iterating through the range of pH from 0-14 in increments of 0.1
      ph = x/100
      charge = self._charge_(ph)
      if abs(charge) < abs(lowest_charge): # comparing calculated pH value and incremented Ph value to find the one closest to zero
        lowest_charge = charge
        lowest_ph = ph
        return lowest_ph # printing out the desired pH

  def aaComposition (self) :
    '''returns the items from dictionary from _initi_ method with keys as single letter amino acid code and the value as the specific count ofeach key'''
    return self.aaComp 

  def _charge_ (self, pH):
    '''calculates the net charge on the protein sequence at a specific pH (a parameter) and the presence of positive and negative charges from particular amino acids in the given
    dictionary'''
    pos_counts =((10 ** ProteinParam.aaNterm)/((10 ** ProteinParam.aaNterm)+(10 ** pH))) # calculating the charge of the N Terminus first and storing it in pos_counts variable
        
    for aa in ProteinParam.aa2chargePos:
      Naa = self.aaComp.get(aa,0)
      pos_charge = ((10 ** ProteinParam.aa2chargePos[aa])/((10 ** ProteinParam.aa2chargePos[aa]) + (10 ** pH)))
      total_pos = pos_charge * Naa # multiplying the charge by the total amount of specific amino acid in the protein sequence
      pos_counts += total_pos # adding each calculated positive charge of specific amino acid into pos_count
            
      neg_counts = ((10 ** pH) /((10 ** ProteinParam.aaCterm) + (10 ** pH))) # calcuating the charge of the C Terminus first and storing it in neg_counts variable
        
    for aa in ProteinParam.aa2chargeNeg: 
      Naa = self.aaComp.get(aa,0)
      neg_charge = ((10 ** pH)/((10 ** ProteinParam.aa2chargeNeg[aa]) + ( 10 ** pH)))
      total_neg = neg_charge * Naa # multiplying the charge by the total amount of specific amino acid in the protein sequence
      neg_counts += total_neg # adding each calculated negative charge of specific amino acid into pos_count
            
      protein_charge = pos_counts - neg_counts # getting the final protein charge by subtracting total negative charges from total positive charges
        
      return protein_charge
                
  def molarExtinction (self,Cystine = True): # Extra Credit for inclusion and exclusion of Cystine
    '''indicates how much light a protein absorbs at a certain wavelength with the presence of specific amino acids within the protein sequence'''
    N_Y = self.aaComp.get("Y",0) # Finding number of specific amino acid within the protein sequence
    N_W = self.aaComp.get("W",0)
    N_C = self.aaComp.get("C",0)
    E_Y = ProteinParam.aa2abs280["Y"] # Extracting the value from the key (amino acid)
    E_W = ProteinParam.aa2abs280["W"]
    E_C = ProteinParam.aa2abs280["C"]
    E = (N_Y*E_Y) + (N_W*E_W) + (N_C*E_C)
    if Cystine:
      return E
    else:
      no_cE = E - (N_C*E_C)
      return no_cE


  def massExtinction (self,Cystine = True): # Extra Credit for inclusion and exclusion of Cystine
    '''Mass extinction coefficient is calculated by the calculated molar extinction divided by the molecular weight of the given protein sequence'''
    if Cystine:
      myMW =  self.molecularWeight()
      return self.molarExtinction() / myMW if myMW else 0.0
    else:
      no_cmw = self.molecularWeight() - (ProteinParam.aa2mw["C"] - ProteinParam.mwH2O)
      return self.molarExtinction(Cystine=False) / no_cmw if no_cmw else 0.0

  def molecularWeight (self):
    '''calculates the molecular weight of the protein sequence and excluding the molecular weight of water formed between the bond of two amino acids'''
    mwcounts = []
    if self.aaCount() > 0:
      for aa in self.protein:
        aa_weight = (ProteinParam.aa2mw[aa] - ProteinParam.mwH2O)
        mwcounts.append(aa_weight) # adding all the calculated mw for each amino acid into the mwcounts list
        total_weight = sum(mwcounts)+ ProteinParam.mwH2O # Adding back the additional molecular weight of water which was subtracted before
        return total_weight
