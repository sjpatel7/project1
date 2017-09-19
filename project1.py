# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 19:26:14 2017

@author: Shawn
"""

import numpy as np

global H1
H1=np.array([[1,1],[1,-1]])/np.sqrt(2)

def H2(theta):
    Q=np.array([[1,0],[0,np.e**(theta*1j)]])
    return Q

global H3  #CNot matrix with cw=ow-1
H3=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])

global H4  #CNot matrix with cw=ow+1
H4=np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]])

global I    #2x2 Identity
I=np.array([[1,0],[0,1]])

global S   #Swapping two wires
S=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])

def HadamardArray(i, k):
# this should apply Hadamard to wire i out of k wires
    n=0 #counter
    if i==0: #if gate on wire 0 
        x=H1
        n=n+1
        while n<k:          #while number of wires tensor multiplied<total number of wires
            x=np.kron(x,I)
            n=n+1
    if i>0:
        x=I
        while n<(i-1):
            x=np.kron(x,I)
            n=n+1
        x=np.kron(x,H1)
        n=n+1
        while n<(k-1):
            x=np.kron(x,I)
            n=n+1

    return x

def PhaseArray(i, k,theta):
# this should apply Phase theta to wire i out of k wires
    n=0 #counter
    if i==0: #if gate on wire 0 
        x=H2(theta)
        n=n+1
        while n<k:          #while number of wires tensor multiplied<total number of wires
            x=np.kron(x,I)
            n=n+1
    if i>0:
        x=I
        while n<(i-1):
            x=np.kron(x,I)
            n=n+1
        x=np.kron(x,H2(theta))
        n=n+1
        while n<(k-1):
            x=np.kron(x,I)
            n=n+1
    

    return x

def CNOTArray(controlWire,otherWire,totalWires):
    cw=controlWire
    ow=otherWire
    tw=totalWires-1
    
    #initial wire number counter
    n=0

    y=0 #initial condition if no swaps or operations have been made
#if control wire is more than one above the other wire
    if cw<ow-1:
        if cw==0:
            x=S
            n=n+1
            while n<tw:
                x=np.kron(x,I)
                n=n+1
            y=x
            cw=cw+1
            n=0
        while cw<ow-1:
            x=I
            while n<cw-1:
                x=np.kron(x,I)
                n=n+1
            x=np.kron(x,S)
            n=n+2
            while n<tw:
                x=np.kron(x,I)
                n=n+1
            if np.size(y)==1:
                y=x
            else:
                y=np.dot(x,y)
            cw=cw+1
            n=0

#if control wire is one above (in diagram) the other wire
    if cw==(ow-1):      
        if cw==0:
            x=H3
            n=n+1
            while n<tw:
                x=np.kron(x,I)
                n=n+1
        if cw>0:
            x=I
            while n<(cw-1):
                x=np.kron(x,I)
                n=n+1
            x=np.kron(x,H3)
            n=n+2
            while n<tw:
                x=np.kron(x,I)
                n=n+1
            if np.size(y)==1:
                y=x
            else:
                y=np.dot(x,y)

#Swap wires back if they were swapped
    while cw!=controlWire:
            x=I
            while n<cw-2:
                x=np.kron(x,I)
                n=n+1
            x=np.kron(x,S)
            n=n+2
            while n<tw:
                x=np.kron(x,I)
                n=n+1
            if np.size(y)==1:
                y=x
            else:
                y=np.dot(x,y)
            cw=cw-1
            n=0

#if control wire is more than one below the other wire
    while cw>ow+1:
        x=I
        while n<cw-2:
            x=np.kron(x,I)
            n=n+1
        x=np.kron(x,S)
        n=n+2
        while n<tw:
            x=np.kron(x,I)
            n=n+1
        if np.size(y)==1:
            y=x
        else:
            y=np.dot(x,y)
        cw=cw-1
        n=0
#if control wire is one below (in diagram) the other wire
    if cw==ow+1:      
        if ow==0:
            x=H4
            n=n+1
            while n<tw:
                x=np.kron(x,I)
                n=n+1
        if ow>0:
            x=I
            while n<(ow-1):
                x=np.kron(x,I)
                n=n+1
            x=np.kron(x,H4)
            n=n+2
            while n<tw:
                x=np.kron(x,I)
                n=n+1
        if np.size(y)==1:
            y=x
        else:
            y=np.dot(x,y)
        n=0
#Swap control wire back if swapped
    while cw!=controlWire:
        x=I
        while n<cw-1:
            x=np.kron(x,I)
            n=n+1
        x=np.kron(x,S)
        n=n+2
        while n<tw:
            x=np.kron(x,I)
            n=n+1
        if np.size(y)==1:
            y=x
        else:
            y=np.dot(x,y)
        cw=cw+1
        n=0
            
    return y


def ReadInput(fileName):
    myInput_lines=open(fileName).readlines()
    myInput=[]
    numberOfWires=int(myInput_lines[0])
    for line in myInput_lines[1:]:
        myInput.append(line.split())
    return (numberOfWires,myInput)
myInput=ReadInput("MyGateDescription.txt")
numWires=myInput[0]
myInput=myInput[1]

i=0
Measure=False
while i<len(myInput):
    if myInput[i][0]=='Measure':
        Measure=True
        break
    gate=myInput[i][0]
    gateParams=np.array([])
    j=1
    while j<len(myInput[i]):     
        gateParams=np.append(gateParams,myInput[i][j])
        j=j+1
    if i==0:
        if gate=='H':
            matrix=HadamardArray(int(gateParams[0]),numWires)
        if gate=='P':
            matrix=PhaseArray(int(gateParams[0]),numWires,float(gateParams[1]))
        if gate=='CNOT':
            matrix=CNOTArray(int(gateParams[0]),int(gateParams[1]),numWires)
        x=matrix
    if i>0:
        if gate=='H':
            matrix=HadamardArray(int(gateParams[0]),numWires)
        if gate=='P':
            matrix=PhaseArray(int(gateParams[0]),numWires,float(gateParams[1]))
        if gate=='CNOT':
            matrix=CNOTArray(int(gateParams[0]),int(gateParams[1]),numWires)
        x=np.dot(matrix,x)
    i=i+1

#function to parse file myInputState.txt
def ReadInput1(fileName):
    myInput_lines=open(fileName).readlines()
    myInput=[]
    for line in myInput_lines[0:]:
        myInput.append(line.split())
    return (myInput)

#get input from myInputState.txt
R=ReadInput1("myInputState.txt")
v=np.array([])

#sigfigs is used for a procedure with multiplying float numbers late
sigfigs=0
for i in range(len(R)):
    v=np.append(v,complex(float(R[i][0]),float(R[i][1])))
    if len(R[i][0])>sigfigs or len(R[i][1])>sigfigs:
        if len(R[i][0])>len(R[i][1]):
            sigfigs=len(R[i][0])-2
        else:
            sigfigs=len(R[i][1])-2
            
v=np.reshape(v,(2**numWires,1))

a=np.dot(x,v)
a=np.multiply(a,10**sigfigs)
l=np.dot(np.reshape(a,(1,8)),np.conjugate(a))
l=np.multiply(l,10**(-2*sigfigs))
l=round(l[0][0],2) 
print("vv*_dagger: ",l)
print()
a=np.multiply(a,np.conj(a))
a=np.multiply(a,10**(-2*sigfigs))
a=a.real
a=a.round(3)



print("output vector:")
print(a)
print()
#randomly measuring
if Measure:
    import matplotlib.pyplot as plt
   
    at=np.transpose(a)
    at=at.tolist()
    at=at[0]

    m=np.random.choice(2**numWires,1000,p=at)
    
    plt.hist(m,bins=8,range=(0,2**numWires-1))
    plt.title("measurements from circuit")
    plt.show()
    print()


