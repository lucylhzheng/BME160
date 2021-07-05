#!/usr/bin/env python3 
# Name: Lucy Zheng (lzheng20)
# Group Members: Christopher Tran (chlotran)

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.
Replace block of (N)s with the {count} determined as the output

Example: 
input: AaNNNNNNNGTC
output: AA{7}GTC

Any lower case letters are converted to uppercase
'''

class DNAstring (str):
    
    def purify(self):
        ''' Return an upcased version of the string,count Ns using index function and returns new DNA sequence.'''
        data = self.upper() #input strings are all upper case 
        while "N" in data: #while loop to iterate the Ns and begin count process
            x = data.find("N")
            y = data.rfind("N")
            num_n = (y - x) + 1 #finding the number of Ns
            new_string =(data[:x] + "{"+ "{0}".format(num_n)+ "}" + data[y+1:]) 
            return new_string
            
    
def main():
    ''' Get user DNA data and clean up the Ns .'''
    data = input('DNA data?')
    thisDNA = DNAstring (data) #inserting data into DNAstring class
    pureData = thisDNA.purify() 
    print(pureData)

main()
