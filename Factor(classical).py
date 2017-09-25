# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 19:55:45 2017

@author: Shawn
"""
import sympy
import numpy as np
import math 
from fractions import gcd
def FindPeriod(N):
    rEven = True
    factors=np.array([])
    while rEven:
        x = np.random.randint(2, np.sqrt(N))
        new_N = N
        while math.gcd(x,new_N) != 1 and new_N !=2:  #not coprime condition
            newFactor = math.gcd(x,new_N)
            factors = np.append(factors,newFactor)
            new_N = int(new_N/newFactor)
            factors = np.append(factors,new_N)
        if math.gcd(x,N) == 1: #coprime condition
            r = 1
            while (x**r - 1) % N !=0:
                r = r + 1
            if r % 2 == 1:
                rEven = False
                factor1 = math.gcd(N,int((x**(r/2)-1)%N))
                factor2 = math.gcd(N,int((x**(r/2)+1)%N))
                factors1 = [factor1,factor2]
                factors = np.append(factors,factors1)
        else:
            rEven = False
    return factors

def Factor(N):
    factors = np.array([])
    if sympy.isprime(N) == True:
        print(N," is a prime number")
        factors = np.append(factors,1)
        factors = np.append(factors,N)
    else:
        print("Factors of ",N,":")
        if N % 2 == 0:
            factors = np.append(factors,2)
            evenFactor = N/2
            factors = np.append(factors,evenFactor)
            while (evenFactor % 2) == 0 and evenFactor!=2:
                factors = np.append(factors,evenFactor/2)
                evenFactor = evenFactor/2
                factors = np.append(factors,N/evenFactor)
        i = 2
        while i <= math.log2(N):
            x = N ** (1 / i)
            if math.isclose(x % 1,0):
                factors = np.append(factors,int(x))
            i = i + 1
        
        factors = np.append(factors, FindPeriod(N))

                
    factors = set(factors)
    return factors

