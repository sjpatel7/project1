# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 19:26:14 2017

@author: Shawn
"""

import numpy as np

global H1
H1 = np.array([[1,1],[1,-1]]) / np.sqrt(2)
 
#H2 is the Phase matrix
def H2(theta):
    Q = np.array([[1,0],[0,np.e**(theta*1j)]])
    return Q

global H3  #CNot matrix with cw=ow-1
H3 = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])

global H4  #CNot matrix with cw=ow+1
H4 = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]])

global I    #2x2 Identity
I = np.array([[1,0],[0,1]])

global S   #Swapping two wires (CNOT1*CNOT2*CNOT1)
S = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])

def CPG(theta):
    CPG = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,np.e**(theta*1j)]])
    return CPG

# this should apply Hadamard gate to wire i out of k wires
def HadamardArray(i, k):
    n = 0 #wire we're on (building the matrix one wire at a time)
    if i == 0: #if gate on wire 0 
        matrix = H1  
        n = n + 1    #wire is updated to the next wire (which hasn't been included in the matrix yet)
        while n < k:          
            matrix = np.kron(matrix,I)
            n = n + 1
    if i > 0:
        matrix = I
        n = n + 1
        while n < i:    
            matrix = np.kron(matrix,I)
            n = n + 1
        matrix = np.kron(matrix,H1)
        n = n + 1
        while n < k:
            matrix = np.kron(matrix,I)
            n = n + 1
    return matrix

#Inverse gates for checkpoint 4 (random circuit)
def InvHadamardArray(i,k):
    invHad = np.linalg.inv(HadamardArray(i,k))
    return invHad
def InvPhaseArray(i,k,theta):
    invPhase = np.linalg.inv(PhaseArray(i,k,theta))
    return invPhase
def InvCNOTArray(cw,ow,tw):
    invCNOT= np.linalg.inv(CNOTArray(cw,ow,tw))
    return invCNOT
def InvCPGate(cw,ow,theta,tw):
    invCPGate = np.linalg.inv(CPGate(cw,ow,theta,tw))
    return invCPGate

# this should apply Phase theta to wire i out of k wires
def PhaseArray(i, k,theta):
    n = 0 #wire we're on 
    if i == 0: #if gate on wire 0 
        matrix = H2(theta)
        n = n + 1     #wire updated to next one (which hasn't been included in the matrix yet)
        while n < k:          
            matrix = np.kron(matrix,I)
            n = n + 1
    if i > 0:
        matrix = I
        n = n + 1
        while n < i:
            matrix = np.kron(matrix,I)
            n = n + 1
        matrix = np.kron(matrix,H2(theta))
        n = n + 1
        while n < k:
            matrix = np.kron(matrix,I)
            n = n + 1
    return matrix

def CNOTArray(controlWire,otherWire,totalWires):
    cw = controlWire
    ow = otherWire
    tw = totalWires
    #wire we are on
    n = 0
    noSwaps = True #initial condition if no swaps or operations have been made
    matrix = np.eye(2**tw)  #this is the resulting matrix of the CNOT, initialized to identity

    if cw < ow - 1:
        matrix = Swap(cw,ow - 1,tw)
        cw = ow - 1
        noSwaps = False
    
    if cw > ow + 1:
        matrix = Swap(cw,ow + 1,tw)
        cw = ow + 1
        noSwaps = False
             
#if control wire is one above (in diagram) the other wire
    matrix_RN = np.eye(2**tw)
    if cw == ow + 1:
        if ow == 0:
            matrix_RN = H4
            n = n + 2
            while n < tw:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
        if ow > 0:
            matrix_RN = I
            n = n + 1
            while n < ow:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
            matrix_RN = np.kron(matrix_RN,H4)
            n = n + 2
            while n < tw:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
    if cw == (ow - 1):      
        if cw == 0:
            matrix_RN = H3
            n = n + 2
            while n < tw:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
        if cw > 0:
            matrix_RN = I
            n = n + 1
            while n < cw:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
            matrix_RN = np.kron(matrix_RN,H3)
            n = n + 2
            while n < tw:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
    if noSwaps:
        matrix = matrix_RN
    else:
        matrix = np.dot(matrix,np.dot(matrix_RN,matrix))
    return matrix

def Swap(firstWire,secondWire,totalWires):
    #I am going to make fw be the wire above and sw be the wire below
    if firstWire < secondWire:
        fw = firstWire
        sw = secondWire
    else:
        sw = firstWire
        fw = secondWire 
        
    fw1 = fw  #initial value of the higher wire
    tw = totalWires
    #initial wire number counter
    n = 0
    matrix = np.eye(2**tw)
    if fw == 0:
        matrix_RN = S
        n = n + 2
        while n < tw:
            matrix_RN = np.kron(matrix_RN,I)
            n = n + 1
        fw = fw + 1
        matrix = matrix_RN
        n = 0
    while fw < sw:
        matrix_RN = I
        n = n + 1
        while n < fw:
            matrix_RN = np.kron(matrix_RN,I)
            n = n + 1
        matrix_RN = np.kron(matrix_RN,S)
        n = n + 2
        while n < tw:
            matrix_RN = np.kron(matrix_RN,I)
            n = n + 1
        matrix = np.dot(matrix_RN,matrix)
        fw = fw + 1
        n = 0
    sw = sw - 1
    
    if sw > 0: #if sw is zero, then the whole swap is complete (wires 0 and 1 were switched)
        while sw != fw1:
            n = 0
            matrix_RN = I
            n = n + 1
            while n < sw - 1:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
            matrix_RN = np.kron(matrix_RN,S)
            n = n + 2
            while n < tw:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
            matrix = np.dot(matrix_RN,matrix)
            sw = sw - 1
    return matrix

#controlled-phase gate
def CPGate(cw,ow,theta,tw):
    #wire we are on
    n = 0
    noSwaps = True #initial condition if no swaps or operations have been made
    matrix = np.eye(2**tw)  #this is the resulting matrix of the CNOT, initialized to identity

    if cw < ow - 1:
        matrix = Swap(cw,ow - 1,tw)
        cw = ow - 1
        noSwaps = False
    
    if cw > ow + 1:
        matrix = Swap(cw,ow + 1,tw)
        cw = ow + 1
        noSwaps = False
             
#if control wire is one above (in diagram) the other wire
    matrix_RN = np.eye(2**tw)
    if cw == ow + 1:
        if ow == 0:
            matrix_RN = CPG(theta)
            n = n + 2
            while n < tw:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
        if ow > 0:
            matrix_RN = I
            n = n + 1
            while n < ow:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
            matrix_RN = np.kron(matrix_RN,CPG(theta))
            n = n + 2
            while n < tw:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
    if cw == (ow - 1):      
        if cw == 0:
            matrix_RN = CPG(theta)
            n = n + 2
            while n < tw:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
        if cw > 0:
            matrix_RN = I
            n = n + 1
            while n < cw:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
            matrix_RN = np.kron(matrix_RN,CPG(theta))
            n = n + 2
            while n < tw:
                matrix_RN = np.kron(matrix_RN,I)
                n = n + 1
    if noSwaps:
        matrix = matrix_RN
    else:
        matrix = np.dot(matrix,np.dot(matrix_RN,matrix))
    return matrix
    
#random ciruit generator (this was used for testing various circuits)
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
#RandomCircuit()
x=CNOTArray(2,0,3)

def ReadInput(fileName):
    myInput_lines = open(fileName).readlines()
    myInput = []
    numberOfWires = int(myInput_lines[0])
    for line in myInput_lines[1:]:
        myInput.append(line.split())
    return (numberOfWires,myInput)
myInput = ReadInput("MyGateDescription.txt")
numWires = myInput[0]
myInput = myInput[1]

i = 0
Measure = False
#this gate builds the matrix from the circuit file 
while i < len(myInput):
    if myInput[i][0] == 'Measure':
        if i == 0: #if no gates
            x = np.eye(2**numWires)            
        Measure = True
        break
    gate = myInput[i][0]  #first character on line specifies gate
    gateParams = np.array([])
    j = 1
    while j < len(myInput[i]):     
        gateParams = np.append(gateParams,myInput[i][j])
        j = j + 1
    if i == 0: #if this is the first gate
        if gate == 'H':
            matrix = HadamardArray(int(gateParams[0]),numWires)
        if gate == 'P':
            matrix = PhaseArray(int(gateParams[0]),numWires,float(gateParams[1]))
        if gate == 'CNOT':
            matrix = CNOTArray(int(gateParams[0]),int(gateParams[1]),numWires)
        if gate == 'CPG':
            matrix = CPGate(int(gateParams[0]),int(gateParams[1]),float(gateParams[2]),numWires)
        x=matrix
    if i>0:
        if gate=='H':
            matrix=HadamardArray(int(gateParams[0]),numWires)
        if gate=='P':
            matrix=PhaseArray(int(gateParams[0]),numWires,float(gateParams[1]))
        if gate=='CNOT':
            matrix=CNOTArray(int(gateParams[0]),int(gateParams[1]),numWires)
        if gate == 'CPG':
            matrix = CPGate(int(gateParams[0]),int(gateParams[1]),float(gateParams[2]),numWires)
            
        if gate=='H_i':
            matrix = InvHadamardArray(int(gateParams[0]),numWires)
        if gate=='P_i':
            matrix = InvPhaseArray(int(gateParams[0]),numWires,float(gateParams[1]))
        if gate == 'CNOT_i':
            matrix= InvCNOTArray(int(gateParams[0]),int(gateParams[1]),numWires)
        if gate == 'CPG_i':
            matrix = InvCPGate(int(gateParams[0]),int(gateParams[1]),float(gateParams[2]),numWires)
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

Swap(0,2,3)
a=np.dot(x,v)
#print("output of circuit with inverse circuit added\n")
#print(a.round(5))
#print()
a=np.multiply(a,10**sigfigs)
l=np.dot(np.reshape(a,(1,8)),np.conjugate(a))
l=np.multiply(l,10**(-2*sigfigs))
l=round(l[0][0],3) 
print("vv*_dagger: ",l)
print()
a=np.multiply(a,np.conj(a))
a=np.multiply(a,10**(-2*sigfigs))
a=a.real
a=a.round(5)

at=np.transpose(a)
print("output vector:")
print(a)
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

#QFT
myInput = ReadInput("QFT.txt")
numWires = myInput[0]
myInput = myInput[1]

i = 0
Measure = False
#this gate builds the matrix from the circuit file 
while i < len(myInput):
    gate = myInput[i][0]  #first character on line specifies gate
    gateParams = np.array([])
    j = 1
    while j < len(myInput[i]):     
        gateParams = np.append(gateParams,myInput[i][j])
        j = j + 1
    if i == 0: #if this is the first gate
        if gate == 'H':
            matrix = HadamardArray(int(gateParams[0]),numWires)
        if gate == 'CPG':
            matrix = CPGate(int(gateParams[0]),int(gateParams[1]),float(gateParams[2]),numWires)
        x1 = matrix
    if i>0:
        if gate=='H':
            matrix=HadamardArray(int(gateParams[0]),numWires)
        if gate == 'CPG':
            matrix = CPGate(int(gateParams[0]),int(gateParams[1]),float(gateParams[2]),numWires)
        x1 = np.dot(matrix,x1)
    i=i+1

output = np.dot(x1,v)
print("Output of QFT")
print(output)