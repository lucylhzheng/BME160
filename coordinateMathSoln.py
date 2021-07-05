#!/usr/bin/env python3 
# Name: Lucy Zheng (lzheng20)


'''
Three sets of atomic coordinates are taken to generate respective bond lengths and bond angles using programmed mathemtical functions

example:
input: C = (39.447, 94.657, 11.824) N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465)
output: 
N-C bond length = 1.33
N-Ca bond length = 1.46
C-N-Ca bond angle = 124.0

make sure the bond angle is outputted as degrees
'''

import math
class Triad :
    """
    Calculate angles and distances among a triad of points.
 
    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Construct a Triad. 
        
        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ). 
            make sure that the numbers within p, q, and r are tuples
        """
        self.p = p
        self.q = q
        self.r = r  
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) /   math.sqrt(self.d2(self.q,self.p)*self.d2(self.r,self.p)))
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /  math.sqrt(self.d2(self.p,self.q)*self.d2(self.r,self.q)))
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /  math.sqrt(self.d2(self.p,self.r)*self.d2(self.q,self.r)))

def main():
    '''Atomic coordinates are collected and separated into each specific tuples after float conversion as it goes through mathematical calculations
    from Triad class'''
    coords = input("Enter three sets of coordinates:")
    replace_leftpara = coords.replace("(","") #replaces parenthesis with commas
    replace_rightpara = replace_leftpara.replace(")",",")
    c1 = float(replace_rightpara[4:10])
    c2 = float(replace_rightpara[12:18])
    c3 = float(replace_rightpara[20:26])
    c = (c1,c2,c3) #tuple of carbon atom
    n1 = float(replace_rightpara[32:38])
    n2 = float(replace_rightpara[40:46])
    n3 = float(replace_rightpara[48:54])
    n = (n1,n2,n3) #tuple of Nitrogen atom
    ca1 = float(replace_rightpara[61:67])
    ca2 = float(replace_rightpara[69:75])
    ca3 = float(replace_rightpara[77:83])
    ca = (ca1,ca2,ca3) #tuple of Calcium atom
    points = Triad(c,n,ca) #inserting each tuple into Triad class
    distanceN_C = points.dPQ()
    print("N-C bond length = "'%s'%float('%.3g'%distanceN_C))
    distanceN_Ca = points.dQR()
    print("N-Ca bond length = "'%s'%float('%.3g'% distanceN_Ca))
    angleq = points.angleQ() 
    angle = angleq
    angle_degrees = angle * 180/math.pi #turning radian into degree
    print("C-N-Ca bond angle = "'%s'%float('%.3g'% angle_degrees))
    
main()
