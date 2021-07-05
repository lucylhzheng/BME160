#!/usr/bin/env python3 
# Name: Lucy Zheng (lzheng20)


'''
Seqname of FASTQ file is inserted and each section of the information is parsed to generate
individual fields of outputs

Example:
input: @EAS139:136:FC706VJ:2:2104:15343:197393
output: 
Instrument = EAS139
Run ID = 136
Flow Cell ID = FC706VJ
Flow Cell Lane = 2
Tile Number = 2104
X-coord = 15343
Y-coord = 197393

make sure "@" is not included in the output
'''

class FastqString:
    ''' Each individual information is separated into parameters which will be stored in the init method and output each specific lines in the 
    parse method'''
    def __init__(self,instrument,runid,flow_id,flow_lane,tile_num,x_coor,y_coor):
        self.instrument = instrument
        self.runid = runid
        self.flow_id = flow_id
        self.flow_lane = flow_lane
        self.tile_num = tile_num
        self.x_coor = x_coor
        self.y_coor = y_coor
    def parse(self):
        '''Each field of information is outputted within their individual output'''
        print ("Instrument = {0}".format(self.instrument))
        print ("Run ID = {0}".format(self.runid))
        print ("Flow Cell ID = {0}".format(self.flow_id))
        print ("Flow Cell Lane = {0}".format(self.flow_lane))
        print ("Tile Number = {0}".format(self.tile_num))
        print ("X-coord = {0}".format(self.x_coor))
        print ("Y-coord = {0}".format(self.y_coor))
    
def main():
    '''the data of seqname is collected and separated into individual variables as it enters the FastQString class'''
    data = input("Enter Seqname line:")
    remove_0 = data[1:] #'@' is removed from the rest of the string
    x = remove_0.split(":") #string is split into individual sections
    instrument = x[0] #each section is categorized within a specific variable
    runid = x[1]
    flow_id = x[2]
    flow_lane = x[3]
    tile_num = x[4]
    x_coor = x[5]
    y_coor = x[6]
    sections = FastqString(instrument,runid,flow_id,flow_lane,tile_num,x_coor,y_coor) #each variable goes through the FastQString class
    sections.parse() #calling the parse method
    
    
    

main()
