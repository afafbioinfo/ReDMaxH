import os, subprocess
#from itertools import chain, combinations
import math

def combinations(data):
	cumul=[]
	combine=[]
	Lisa=[elem for elem in data]
	for elem in Lisa:
		cumul.append(elem)
		X= list(cumul) #list to be able to keep the different combinations
		#print X
		combine.append(X)
	return combine
	
	

def overlppingPairhelices(H1,H2):
	#print H1,H2
	s1=H1[1][0]
	e1=H1[1][1]
	s2=H2[1][0]
	e2=H2[1][1]
	L1=H1[1][2]
	L2=H2[1][2]
	if (s2>=e1 ) or (e2<=s1) or (s2>=s1+L1 and e2<=e1-L1) or (s1>=s2+L2 and e1<=e2-L2):
		return False
	else:
		return True
		
def CascadingConstraints(ListHelices):
	n=len(ListHelices)
	Combinations=[]
	#initialize  with the first helix
	C=dict()
	L=0
	C[L]=[ListHelices[0]]
	for k in range(1,n):
		
		#print k ,"gggg"
		for index in C.keys():
			
			#print "index", index
			Lii=[]
			for helix in C[index]:
				
				#print "check",helix,ListHelices[k],overlppingPairhelices(helix,ListHelices[k])
				if overlppingPairhelices(helix,ListHelices[k]):
					Lii.append(helix)
			if len(Lii)==0:#no overlapping
				C[index].append(ListHelices[k])
				
			else:
				#L+=1
				
				C[index]=[ elem for elem in C[index] if elem not in Lii]# rmove the overlapping helices
				C[index].append(ListHelices[k])
			#print "C", C
	for elem in C:
		if C[elem] not in Combinations:
			Combinations.append(C[elem])
	return  Combinations


def CascadingConstraintsversion2(ListHelices):
	n=len(ListHelices)
	HH=[]
	Stack=[[ListHelices[0]]]
	#initialize  with the first helix
	C=dict()
	L=0
	warning=[]
	C[L]=[ListHelices[0]]
	for k in range(1,n):
		
		#print k ,"gggg"
		for index in C.keys():
		
			#print "index", index
			Lii=[]
			for helix in C[index]:
				
				#print "check",helix,ListHelices[k],overlppingPairhelices(helix,ListHelices[k])
				if overlppingPairhelices(helix,ListHelices[k]):
					Lii.append(helix)
			#Stack.append(C[index])
			if len(Lii)==0:#no overlapping
				
				#print "whyyy",Stack
				C[index].append(ListHelices[k]) 
				#if C[index] not in Stack:
				#	#Stack.append(C[index])
			else:
				L+=1
				#Stack.append([ListHelices[k]]) #print out a helix if there is an overlap with the precedent pack of helices
				C[L]=[ elem for elem in C[index] if elem not in Lii]# rmove the overlapping helices
				warning.append( str(ListHelices[k])+ "is not compatible with"+ str(Lii))
				C[L].append(ListHelices[k])
				#Stack.append([ListHelices[k]])
				#if C[L] not in Stack:
				#	#Stack.append(C[L])
			#print "C", C
	#print "ggggggggggC",C
	HH.append([ListHelices[0]])
	for elem in C:
		for k in range(2,len( C[elem])+1):
			if C[elem][:k] not in HH:
				HH.append( C[elem][:k])	
	#print set(warning)
	#print HH
	return HH
	
def ComputeCombinations(ListHelices):
	Combinations=[]
	n=len(ListHelices)
	for i in range(n):
		#print ListHelices[i]
		#Combinations.append((ListHelices[i]))
		for k in range(i+1,n):
			if not overlppingPairhelices(ListHelices[i],ListHelices[k]):
				Combinations.append([ListHelices[i],ListHelices[k]])
				for h in range(k+1,n):
					if not overlppingPairhelices(ListHelices[k],ListHelices[h]) and not overlppingPairhelices(ListHelices[i],ListHelices[h]):
						Combinations.append([ListHelices[i],ListHelices[k],ListHelices[h]])
				
						for M in range(h+1,n):
							if not overlppingPairhelices(ListHelices[h],ListHelices[M]) and not overlppingPairhelices(ListHelices[k],ListHelices[M]) and not overlppingPairhelices(ListHelices[i],ListHelices[M]):
								Combinations.append([ListHelices[i],ListHelices[k],ListHelices[h],ListHelices[M]])
			
	return  Combinations
	
# Parse a file by returning lines it contains
def Parsefile(Path):
    fileIn = open(Path, "r")
    lines = [l.strip() for l in fileIn.readlines()]
    fileIn.close()
    return lines
    
def GetlinefromFile(Path, Linenumber):
    return Parsefile(Path)[Linenumber]
    
def removefileswithextension(dir):
	files_in_directory = os.listdir(dir)
	filtered_files = [file for file in files_in_directory if file.endswith(".ps")]
	for file in filtered_files:
		path_to_file = os.path.join(dir, file)
		os.remove(path_to_file)
    
def buildConstraint(sequence, listofhelices,Method):
	#Method=1: < > as constraints (specific BPs), Method=2: | | as constraints ncts are paired, Method=3 constrain the middle pairing by: removing 1 BP from each extremity if length <7 remove 2bps from extremities otherwise 
	#print "sososooso",listofhelices
	structure=["." for i in sequence]
	for helix in listofhelices:
		#print "dddd",helix
		if 1==1:
			#print helix
			if Method==1:
				for j in range(helix[1][2]):
					#print "op", helix[1][0]+j,"clo", helix[1][1]-j
					structure[helix[1][0]+j-1]="("
					structure[helix[1][1]-j-1]=")"
			if Method==2:
				for j in range(helix[1][2]):
					#print "op", helix[1][0]+j,"clo", helix[1][1]-j
					structure[helix[1][0]+j-1]="|"
					structure[helix[1][1]-j-1]="|"
			if Method==3:
				if helix[1][2]<7:
					for j in range(1,helix[1][2]-1):#one BP from each extremitie
						#print "op", helix[1][0]+j,"clo", helix[1][1]-j
						structure[helix[1][0]+j-1]="("
						structure[helix[1][1]-j-1]=")"
				
				else:
					for j in range(2,helix[1][2]-2):#one BP from each extremitie
						#print "op", helix[1][0]+j,"clo", helix[1][1]-j
						structure[helix[1][0]+j-1]="("
						structure[helix[1][1]-j-1]=")"
						
			if Method==4:
				# Constraint print outprint helix, helix[1][2]
				j=int(math.floor(helix[1][2]*0.5)) #math.floor() function rounds down to the next full integer.
				#print j
				structure[helix[1][0]+j-1]="("
				structure[helix[1][1]-j-1]=")"
			# Constraint print out print helix,"\n",("".join(structure))
			if Method==5:#3 base pairs
				# Constraint print outprint helix, helix[1][2]
				j=int(math.floor(helix[1][2]*0.5)) #math.floor() function rounds down to the next full integer.
				#print j
				structure[helix[1][0]+j-1]="("
				structure[helix[1][1]-j-1]=")"
				structure[helix[1][0]+j-2]="("
				structure[helix[1][1]-j-2]=")"
				structure[helix[1][0]+j]="("
				structure[helix[1][1]-j]=")"
			# Constraint print out print helix,"\n",("".join(structure))
	structure=("".join(structure)) 
	
	outputfile="contraint.txt"
    	o = open(outputfile, "w")  # geneate intermediate file with sequence+constraint
    	o.write(">%s%s\n" % ([helix[0] for helix in listofhelices],"method"+str(Method)))
    	o.write("%s\n%s\n" % (sequence, structure))
    	#print listofhelices,sequence, structure
    	o.close()
    	return outputfile
def Predict2D(Constraintfile,Tag):
        #print Constraintfile   
        Input=Constraintfile
	output="result.txt"
	Command = 'RNAfold -p --noLP -d2 '
	if Tag!="MFE":
      		Command += ' -C --enforceConstraint '
  	#print Command 
	#print "before command line"
	subprocess.call(Command, stdin=open(Input, 'r'), stdout=open(output, 'wb'),
                            stderr=open(os.devnull, 'w'), shell=True)
     	#print GetlinefromFile(output, 0),GetlinefromFile(output, 2)
     	filesize = os.path.getsize(output)
     	#print "lllllllllllllll",filesize,GetlinefromFile(output,2)
	if filesize != 0:
		if Tag!="MFE":
			return (GetlinefromFile(output, 0),GetlinefromFile(Constraintfile, 2),GetlinefromFile(output, 2))
		else:
			#print "strcuture",GetlinefromFile(output, 2)
			return (GetlinefromFile(output, 0),"No constraints",GetlinefromFile(output, 2))
    	else:
    		return "No generated structure for"+ Tag
    		
## assessing predicted structures pseudoknots should be represented by the symbol <<..>>
def ListBasePairsFromStruct(Struct):  # return dic={structure:[liste de pairs de base ],....}
    lista = []
    stack = []
    stackpseudoknots= []
    for i in range(len(Struct)):  # sequence length
        #print Struct
        if Struct[i] == '(':  # Opening base pair
            stack.append(i)
        elif Struct[i] == '<':  # Opening base pair
            stackpseudoknots.append(i)
        elif Struct[i] == '>':  # Closing base pair
            k = stackpseudoknots.pop()
            lista.append((k, i))
        elif Struct[i] == ')':  # Closing base pair
            k = stack.pop()
            lista.append((k, i))
    
    return lista


def TP_Rsample(Correct, Predicted):
        count=0
	for (i,j) in Predicted:
		if (i,j) in Correct or (i,j-1) in Correct or (i,j+1) in Correct	or (i+1,j) in Correct or (i-1,j) in Correct:
			count+=1
	return count 
	
def PPV_Sensitivity_specificityinverse_F1_MCC(Correct,Predicted,n):
	    X=Correct
	    Y= Predicted
	    N= n*(n-1)/2.
	    TP= (len(X)+len(Y)-len(set(X).symmetric_difference(set(Y)) ))/2
	    FP= len(Y)-TP
	    FN= len(X)-TP
	    TN= N - len(Y)-len(X)+TP
	    if (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)!=0:
		return TP_Rsample(Correct, Predicted)/(float(len(Correct))), TP_Rsample(Correct, Predicted)/float(len(Predicted)),1-TN/float(TN+FP),2*TP/float(2*TP+FP+FN),(TP*TN-FP*FN)/float(math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
	    else:
		return 0,0,0,0,0
		
def PPV_sensitivity_Rsample(Correct,Predicted):
        if len(Correct)!=0:
		#Sensitivity=len(set(Correct).intersection(set(Predicted)))/(float(len(Correct)))
		Sensitivity=TP_Rsample(Correct, Predicted)/(float(len(Correct)))
	else:
		Sensitivity=0
	if len(Predicted)!=0:
		#PPV= len(set(Correct).intersection(set(Predicted)))/float(len(Predicted))
		PPV= TP_Rsample(Correct, Predicted)/float(len(Predicted))
	else:
		PPV=0
	return Sensitivity, PPV 		




