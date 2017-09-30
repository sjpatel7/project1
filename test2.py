# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 01:04:49 2017

@author: Shawn
"""
import numpy as np
def simulator(x):
    def CPhase(cw,ow,numWires,theta,inputState):
        cBit = numWires - cw
        oBit = numWires - ow
        
        element = 0
        basisTran = np.zeros((np.size(inputState),1),dtype=complex)
        if cBit == oBit:
            basisTran = inputState
        else:
            while element < np.size(inputState):
                if inputState[element] == 0:
                    pass
                else:
                    #condition for cBit and oBit to be both 1
                    if (element % 2**cBit >= 2**(cBit - 1)) and (element % 2**oBit >= 2**(oBit - 1)):
                        basisTran[element] += np.e**(1j*theta)*inputState[element]
                    else: 
                        basisTran[element] += inputState[element]
                element += 1
        return basisTran
    def HadamardArr(wire,numWires,inputState):
        #wire specifies which bit the gate acts on 
        n = numWires - wire  #this variable is the nth to last bit (the wire's bit)
       
        #this outer loop loops over the basis of the input
        element = 0
        
        basisTran = np.zeros((np.size(inputState),1),dtype=complex)  #curren basis transform in loop, will add to 
        #total transformation
        while element < np.size(inputState):
            if inputState[element] == 0:
                pass
            else:
                if (element % 2**n < 2**(n - 1)):      #condition for nth to last bit to be 0
                    basisTran[element] += 1/np.sqrt(2)*inputState[element]
                    basisTran[element + 2**(n - 1)] += (1/np.sqrt(2)) * inputState[element]
                else:  #nth last bit is 1
                    basisTran[element] -= 1/np.sqrt(2)*inputState[element]
                    basisTran[element - 2**(n - 1)] += (1/np.sqrt(2)) * inputState[element]   
            element += 1
        outputState = basisTran
        return outputState
    
    def CNOTArr(cw,ow,numWires,inputState):
        cBit = numWires - cw
        oBit = numWires - ow
        
        element = 0
        basisTran = np.zeros((np.size(inputState),1),dtype=complex)
        if cBit == oBit:
            basisTran = inputState
        else:
            while element < np.size(inputState):
                if inputState[element] == 0:
                    pass
                else:
                    #condition for control bit to be 0
                    if (element % 2**cBit < 2**(cBit - 1)):
                        basisTran[element] += inputState[element]
                    else: #control bit is 1
                        if (element % 2**oBit < 2**(oBit - 1)): #if other bit is 0
                            basisTran[element + 2**(oBit - 1)] += inputState[element]
                        else:  #other bit is 1
                            basisTran[element - 2**(oBit - 1)] += inputState[element]
                element += 1
        return basisTran
    
    def PhaseArr(wire,numWires,theta,inputState):
        n = numWires - wire
        element = 0
        basisTran = np.zeros((np.size(inputState),1), dtype = complex)
        while element < np.size(inputState):
            if inputState[element] == 0:
                pass
            else:
                if (element % 2**n < 2**(n - 1)): # if wire bit is 0
                    basisTran[element] += inputState[element]
                else:
                    basisTran[element] += np.e**(1j*theta)*inputState[element]
            element += 1
        return basisTran

    def ReadInput(fileName):
        myInput_lines = open(fileName).readlines()
        myInput = []
        numberOfWires = int(myInput_lines[0])
        for line in myInput_lines[1:]:
            myInput.append(line.split())
        return (numberOfWires,myInput)
    myInput = ReadInput("randomCircuit.txt")
    numWires = myInput[0]
    myInput = myInput[1]
    
    i = 0
    Measure = False
    #this gate builds the matrix from the circuit file 
    while i < len(myInput):
        if myInput[i][0] == 'Measure':  
            Measure = True
            break
        gate = myInput[i][0]  #first character on line specifies gate
        gateParams = np.array([])
        j = 1
        while j < len(myInput[i]):     
            gateParams = np.append(gateParams,myInput[i][j])
            j = j + 1

        if gate == 'H':
            x = HadamardArr(int(gateParams[0]),numWires,x)
        if gate == 'P':
            x = PhaseArr(int(gateParams[0]),numWires,float(gateParams[1]),x)
        if gate == 'CNOT':
            x = CNOTArr(int(gateParams[0]),int(gateParams[1]),numWires,x)
        #if gate == 'CPG':
            #   matrix = CPGate(int(gateParams[0]),int(gateParams[1]),float(gateParams[2]),numWires)
        #x=matrix
        i=i+1
    
    #print("output of circuit with inverse circuit added\n")
    #print(a.round(5))
    #print()
    l = np.dot(np.reshape(x,(1,np.size(x))),np.conjugate(x))
    print("vv*_dagger: ",l)
    print()

    x = np.multiply(x,np.conj(x))
    x = x.real
    x=x.round(5)
    
    at=np.transpose(x)
    print("output vector:")
    print(x)
    print()
    #randomly measuring
    if Measure:
        import matplotlib.pyplot as plt
        at=at.tolist()
        at=at[0]
        numBins = 2 ** numWires
        if np.sum(at) != 1:
            remainder = 1 - np.sum(at)
            if remainder!=0:
                print("not unitary!!!!")
            if remainder < 0:
                at[np.argmax(at)] = at[np.argmax(at)] + remainder
            else:    
                at = np.append(at,remainder)
                numBins = 2 ** numWires + 1
            
        measurements=np.random.choice(numBins,10000,p = at)
        
        plt.hist(measurements,bins = numBins,range = (0,2**numWires-1))
        plt.title("measurements from circuit")
        plt.show()
        print()
    return x
def RandomCircuit():
    myInput=open("randomCircuit.txt","w")
    #first line in circuit (number of wires)
    # if number of wires is random, uncomment this line string=np.random.randint(1,9)
    string = 3
    myInput.write("%d\n"%string)
    
    #write all gates
    gateList = ['H','P','CNOT','Measure']
    while string != 'Measure':
        g = np.random.choice(gateList)
        if g == 'H':
            n = np.random.randint(0,string)
            myInput.write("%s %d\n"%(g,n))
        if g == 'P':
            n = np.random.randint(0,string)
            t = 2*np.pi*np.random.random_sample()
            myInput.write("%s %d %f\n"%(g,n,t))
        if g == 'CNOT':
            n = np.random.randint(0,string)
            z = np.random.randint(0,string)
            while z == n:
                z = np.random.randint(0,string)
            myInput.write("%s %d %d\n"%(g,n,z))
        if g == 'Measure':
            myInput.write("Measure")
            string = 'Measure'
    myInput.close()
RandomCircuit()

def CPhase(cw,ow,numWires,theta,inputState):
    cBit = numWires - cw
    oBit = numWires - ow
    
    element = 0
    basisTran = np.zeros((np.size(inputState),1),dtype=complex)
    if cBit == oBit:
        basisTran = inputState
    else:
        while element < np.size(inputState):
            if inputState[element] == 0:
                pass
            else:
                #condition for cBit and oBit to be both 1
                if ((element % 2**cBit >= 2**(cBit - 1)) and (element % 2**oBit >= 2**(oBit - 1))):
                    basisTran[element] += np.e**(1j*theta)*inputState[element]
                else: 
                    basisTran[element] += inputState[element]
            element += 1
    return basisTran
def HadamardArr(wire,numWires,inputState):
    #wire specifies which bit the gate acts on 
    n = numWires - wire  #this variable is the nth to last bit (the wire's bit)
       
    #this outer loop loops over the basis of the input
    element = 0
      
    basisTran = np.zeros((np.size(inputState),1),dtype=complex)  #curren basis transform in loop, will add to 
    #total transformation
    while element < np.size(inputState):
        if inputState[element] == 0:
            pass
        else:
            if (element % 2**n < 2**(n - 1)):      #condition for nth to last bit to be 0
                basisTran[element] += 1/np.sqrt(2)*inputState[element]
                basisTran[element + 2**(n - 1)] += (1/np.sqrt(2)) * inputState[element]
            else:  #nth last bit is 1
                basisTran[element] -= 1/np.sqrt(2)*inputState[element]
                basisTran[element - 2**(n - 1)] += (1/np.sqrt(2)) * inputState[element]   
        element += 1
    outputState = basisTran
    return outputState
    
#this is a test of the matrix
def QFTmatrix(N):
    size = 2**N
    qftmatrix = np.eye(size, dtype=complex)
    w = np.e**(2*np.pi*1j/size)
    for i in range(size):
        for j in range(size):
            qftmatrix[i][j] = w**(i*j)
    qftmatrix = qftmatrix/np.sqrt(size)
    return qftmatrix

def QFT(N,inputState):
    x = inputState
    i = 0
    while i < N:
        j = 0
        while j < i:
            x = CPhase(i,j,N,np.pi/(2**(i-j)),x)
            print("CPhase(",i,",",j,",",N,",",np.pi/(2**(i-j)),")")
            j += 1
        x = HadamardArr(i,N,x)
        print("HadamardArr(",i,",",N,")")
        i += 1
    print()
    return x
'''
n=0
mismatch=0
N=4
while n<20:
    x=np.array([])
    for i in range(2**N):
        x = np.append(x,np.random.rand())
    x = np.transpose([x])
    y = np.around(QFT(N,x),5)
    z = np.around(np.dot(QFTmatrix(N),x),5)
    if np.any(y != z):
        mismatch +=1
    n+=1
print("number of mismatches:",mismatch)
'''
n=0
mismatch=0
N=7
from project1 import CPGate
while n<100:
    x=np.array([])
    for i in range(2**N):
        x = np.append(x,np.random.rand())
    x = np.transpose([x])
    ow = np.random.randint(0,N)
    cw = np.random.randint(0,N)
    theta=np.random.rand()*np.pi*2
    while ow == cw:
        ow = np.random.randint(0,N)
    y = np.around(CPhase(cw,ow,N,theta,x),5)
    z = np.around(np.dot(CPGate(cw,ow,theta,N),x),5)
    if np.any(y != z):
        mismatch +=1
        print("(cw,ow): m ", cw,ow)

    n+=1
    print("(cw,ow) ",cw,ow)
print("number of mismatches:",mismatch)
'''
    #this is a test of the matrix
    n = 0
    mismatch=0
    from project1 import HadamardArray
    N=4
    while n<20:   
        x = np.array([])
        for i in range(2**N):
            x = np.append(x,np.random.rand())
            x = np.transpose(x)
        wire = np.random.randint(0,N)
    
        a = np.around(HadamardArr(wire, N, x),5)
        b = np.transpose([np.around(np.dot(HadamardArray(wire,N),x),5)])
        if np.any(a!=b):
            print("mismatch")
            mismatch += 1
        n+=1
    print("number of mismatches",mismatch)
        
    n = 0
    mismatch=0
    from project1 import CNOTArray
    N = 4
    while n<20:   
        x = np.array([])
        for i in range(2**N):
            x = np.append(x,np.random.rand())
            x = np.transpose(x)
        cwire = np.random.randint(0,N)
        owire = np.random.randint(0,N)
        c = np.transpose(np.around(CNOTArr(cwire,owire, N, x),5))
        d = np.transpose([np.around(np.dot(CNOTArray(cwire,owire,N),x),5)])
        if np.any(c!=d):
            print("mismatch, (cw,ow)=",cwire,owire)
            print("c\n",c)
            print("d\n",d)
            mismatch += 1
        n+=1
    print("number of mismatches",mismatch)
    
    n = 0
    mismatch=0
    from project1 import PhaseArray
    N = 4
    while n<20:   
        x = np.array([])
        for i in range(2**N):
            x = np.append(x,np.random.rand())
            x = np.transpose(x)
        wire = np.random.randint(0,N)
        theta = 2*np.pi*np.random.rand()
        y = np.transpose([np.around(PhaseArr(wire, N, theta, x),5)])
        z = np.transpose([np.around(np.dot(PhaseArray(wire,N,theta),x),5)])
        if np.any(y!=z):
            print("mismatch")
            print("y\n",y)
            print("z\n",z)
            mismatch += 1
        n+=1
    print("number of mismatches",mismatch)
    '''