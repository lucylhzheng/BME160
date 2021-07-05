# !/usr/bin/env python3
# Name: Lucy Zheng (lzheng20)
# Group Members: Christopher Tran (chlotran)
'''
Takes in a string of amino acids and calculates the physical-chemical properties related to the given protein sequence.


Example:

Input: VLSPADKTNVKAAW

Output:
Number of Amino Acids: 14
Molecular Weight: 1499.7
molar Extinction coefficient: 5500.00
mass Extinction coefficient: 3.67
Theoretical pI: 9.88
Amino acid composition:
A = 21.43%
C = 0.00%
D = 7.14%
E = 0.00%
F = 0.00%
G = 0.00%
H = 0.00%
I = 0.00%
K = 14.29%
L = 7.14%
M = 0.00%
N = 7.14%
P = 7.14%
Q = 0.00%
R = 0.00%
S = 7.14%
T = 7.14%
V = 14.29%
W = 7.14%
Y = 0.00%


'''
class ProteinParam:

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
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
            
            
# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        for key in keys :
            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))
        inString = input('protein sequence?')
        

if __name__ == "__main__":
    main()
