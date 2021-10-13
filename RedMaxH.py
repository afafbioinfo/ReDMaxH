# -*- coding: utf-8 -*-
import numpy as np
import math
from collections import defaultdict
import glob, os, subprocess
from itertools import chain
import argparse
import operator
import csv 
import Prediction as PR
import re
from sklearn.datasets import make_blobs
regex = re.compile('[^a-zA-Z]')
parser = argparse.ArgumentParser()                                               
parser.add_argument("--file", "-f", type=str, required=True)
parser.add_argument("--profile", "-p", type=str, required=True)
parser.add_argument("--output", "-o", type=str, required=True)
parser.add_argument("--fasta", "-R", type=str, required=True)
parser.add_argument("--amplicon", "-a",  type=int, required=True)
parser.add_argument("--backwardamplicon", "-b",  type=int, required=True)
parser.add_argument("--Helixnbr", "-x",  type=int, required=True)
args = parser.parse_args()

def Merge(dict1, dict2):
    return(dict2.update(dict1))
    
def Computecorrelations(MM,UU,Prior,Pairs,lenRNA,Name,MinLenHelix,setofhelices):
	Rd= np.zeros((lenRNA,lenRNA)) 
	for i in range( len(MM)):
		for j in range(len(MM)):
			if Prior[i,j]!=0:
				Rd[i,j]=(MM[i,j]-Prior[i,j])/float(Prior[i,j])
	#set diagonals to Zero
        for i in range( len(MM)):
		MM[i,i]=np.nan
		Rd[i,i]=0
	
	Rd[Rd==0]=np.nan # set null normalized values to NA because Rd was initially defined as matrix of zeros
	######### Output MM count and Rd
	with open (os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+"_Rawdata.csv"), 'w') as outfile: 
		outfile.write("%s,%s,%s,%s,%s\n"% ( "i","j","ncts","MM","Rd" ))
		for i in range( AmpliconF, lenRNA-AmpliconB):
			for j in range( i+7,lenRNA-AmpliconB):
				outfile.write("%i,%i,%s,%.4f,%.4f \n"% (i+1,j+1,str(seq[i]+seq[j]),round(MM[i,j],4),round(Rd[i,j],4)))
	outfile.close()
	######## output data for specific helices
	outputcorrelation=os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+"_Norm_helices.csv")
	outputcorrelation1=os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+"_RedmaxH_Rankedhelices.csv")
	contributingpairs=os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+"_contributingpairsforrankedhelices.csv")
	Lista=[]
	#print "pairs in helices with length L\n" 
	for Helix in setofhelices:
		if Helix[1][2]>MinLenHelix: 	
			listtuple=[]	
			i=Helix[1][0]-1
			j=Helix[1][1]-1
			l=Helix[1][2]
			###### fill up listtuple with tuples along with their Rd value
			for k in range(l):# neigborhood of order 1 [i-1,i+1]x[j-1,j+1]
				for (n,t) in [(i+k,j-k),(i+k,j-k-1),(i+k,j-k+1),(i+k-1,j-k),(i+k-1,j-k-1),(i+k-1,j-k+1),(i+k+1,j-k),(i+k+1,j-k-1),(i+k+1,j-k+1)]:
					if j+1<lenRNA and math.isnan(Rd[n,t])!=True and (n,t,Rd[n,t])not in listtuple:  
						Lista.append((Helix[0],n,t,Rd[n,t],l))
						listtuple.append((n,t,Rd[n,t]))		
        # Get the #of common pairs between helices NA not included.
	dicto= defaultdict()
	dictobis= defaultdict() #one additional element that is the length of the helix
	CommonPairs=defaultdict()
	Matrixcommon=defaultdict(lambda: defaultdict(lambda:0))
	#initialization
	for elem in Lista:
		dicto[elem[0]]=list()
		dictobis[elem[0]]=list()
		CommonPairs[elem[0]]=list()	
	#incrementation
	for elem in Lista:
		dicto[elem[0]].append((elem[1],elem[2],elem[3]))
		dictobis[elem[0]].append((elem[1],elem[2],elem[3],elem[4]))
	Averagehelix=defaultdict() 
	Averagehelixtoadd=defaultdict() 
	dict3=defaultdict()
	for key in dictobis:
		toaverage=[]
		## CC content and AA content 
		listAA=0
		listCC=0
		for elem in dictobis[key]:
			if seq[elem[0]].upper()=="C" and seq[elem[1]].upper()=="C":
				listCC+=1
			if seq[elem[0]].upper()=="A" and seq[elem[1]].upper()=="A":
				listAA+=1 

			#### Resize the interval to [-1,3]
			if elem[2]>3:
				toaverage.append(3)
			else:
				toaverage.append(elem[2])			
		#tuple of average, nbr of observations, expected number of observations
		Averagehelix[key]=(round(np.mean(toaverage),3),len(toaverage),dictobis[key][0][3]*5+4,listCC,listAA,toaverage,round(np.std(toaverage),3),round(np.median(toaverage),3),round(np.mean([x for x in toaverage if x>0]),3),round(np.std([x for x in toaverage if x>0]),3),round(np.median([x for x in toaverage if x>0]),3))
	#2 ------------ Combine helices = forming helix compounds
	for key1 in Averagehelix:
		firstlist=[elem for elem in dicto[key1]]
		for key2 in Averagehelix:
			if key2!=key1:
				secondlist=[elem for elem in dicto[key2]]
				if len(set(firstlist) & set(secondlist))!=0:
					CommonPairs[key1].append((key2,len(set(firstlist) & set(secondlist))))
				Matrixcommon[key1][key2]=len(set(firstlist) & set(secondlist))	
	removed=[]
	Commonpairsfile="2DHelicesCommonPairs.csv"
	with open (os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+Commonpairsfile), 'w') as outfile:
		outfile.write("%s\n"% ( str([Helix[0] for Helix in setofhelices])))#,"commonpairs" ))
		for i in  [Helix[0] for Helix in setofhelices]:
			outfile.write("%s"% (i))
			for j in [Helix[0] for Helix in setofhelices]:
				outfile.write(",%i"% (Matrixcommon[str(i)][str(j)]))
			outfile.write("\n")
		outfile.close()
	Commonpairsfile="MergedHelices.csv"
	with open (os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+Commonpairsfile), 'w') as outfile:	
	        outfile.write("%s\t%s\t%s\t%s\n"% ("Helix 1","Helix 2","# common pairs","Rc"))
		#combine 
		Donehelices=[]
		for Helix in Averagehelix:
			for Helixrelated in CommonPairs[ Helix]:
				X=[]
				Helice=Helixrelated[0]	
				commonpairs=Helixrelated[1]
				X=Averagehelix[Helix][5]
				X= Averagehelix[Helix][5]+Averagehelix[Helice][5]
				for elem in list(set([elem for elem in dicto[Helix]]) & set([elem for elem in dicto[Helice]])):
					if elem[2]>3:
						X.remove(3) #romove one instance of a value if greater than 3 
					else:
						X.remove(elem[2]) # remove will romove one instance of a value
				if (Helix,Helice) not in Donehelices:
					Donehelices.append((Helix,Helice))
					Donehelices.append((Helice,Helix))
					if  commonpairs>9:#minimum 10 tuples
						outfile.write("%s\t%s\t%i\t%.3f\n"% (Helix,Helice,commonpairs,round(np.mean(X),3)))
						Averagehelixtoadd[str(Helix+"-"+Helice)]=(round(np.mean(X),3),len(X),999,999,999,X,round(np.std(X),3),round(np.median(toaverage),3),round(np.mean([x for x in toaverage if x>0]),3),round(np.std([x for x in toaverage if x>0]),3),round(np.median([x for x in toaverage if x>0]),3))
	
						#feed the list of Rc with the new combinaisons and remove the preexisting ones.
						removed.append(Helix)
						removed.append(Helice)
	outfile.close()
	for Helix in list(set(removed)): # To remove redundant helices
		if Averagehelix[Helix]:
			del Averagehelix[Helix]#
	# remove Underrepresented helices 
	Underrepresentedhelices=[]	
	for helix in Averagehelix:
		if Averagehelix[helix][1]<(3*5+4):
			Underrepresentedhelices.append(helix)
	for helix in list(set(	Underrepresentedhelices)):
		del Averagehelix[helix]

	# combine two sets
	for k, v in chain(Averagehelix.items(),Averagehelixtoadd.items()):
    		dict3[k]=v
	#sort the dictionary:
	SortedAverageHelix = dict3.items()
	SortedAverageHelix.sort(key=lambda x:-x[1][0])			
	with open (outputcorrelation, 'w') as outfile:
		outfile.write("%s,%s,%s,%s,%s,%s,%s\n"% ( "Dataset","Helix","i","j","seqij","Rltdiff","ijdist"))#,"commonpairs" ))
		for elem in Lista:
			outfile.write("%s,%s,%s,%s,%s,%.4f,%i \n"% ( ((Name.split("/")[-1]).split("_")[1]).split(".")[0],elem[0],str(elem[1]+1),str(elem[2]+1),str(seq[elem[1]]+seq[elem[2]]),round(elem[3],4),elem[2]-elem[1]+1))#,CommonPairs[elem[0]]))
	
	Order=[]
	with open (outputcorrelation1, 'w') as outfile:
		#here the number of clusters
		outfile.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"% ( "Helix","averageRd","nbrobservations","Length","nbrCCpairs","nbrAApairs","std","median","averageonpositive","stdonpositive","medianonpsitive" ))
		for elem in SortedAverageHelix:
			Order.append(elem[0])
			outfile.write("%s,%.4f,%i,%i,%i,%i,%.4f,%.4f,%.4f,%.4f,%.4f\n"% (elem[0], elem[1][0] , elem[1][1], (elem[1][2]-4)/5, elem[1][3], elem[1][4], elem[1][6], elem[1][7], elem[1][8], elem[1][9], elem[1][10] ))
			#get profiling
			#if '-' not in elem[0]:
			#	outfile.write("%s,%.4f,%i,%i,%i,%i,%s,%.4f,%.4f,%.4f,%.4f,%.4f\n"% (elem[0], elem[1][0] , elem[1][1], (elem[1][2]-4)/5, elem[1][3], elem[1][4], [i[2] for i in setofhelices if i[0]==elem[0]][0], elem[1][6], elem[1][7], elem[1][8], elem[1][9], elem[1][10] ))
			#else:
			# a combined profile
			#	outfile.write("%s,%.4f,%i,%i,%i,%i,%s,%.4f,%.4f,%.4f,%.4f,%.4f\n"% (elem[0], elem[1][0] , elem[1][1], (elem[1][2]-4)/5, elem[1][3], elem[1][4],  [i[2] for i in setofhelices if i[0]==elem[0].split("-")[0]][0] + " &  "+ [i[2] for i in setofhelices if i[0]==elem[0].split("-")[1]][0] , elem[1][6], elem[1][7], elem[1][8], elem[1][9], elem[1][10]))
	outfile.close()
	with open(contributingpairs,'w') as outfile:
		outfile.write("%s,%s\n"%("Helix","Rforpairs"))
		for elem in SortedAverageHelix:
				#print elem[1][5]
				for k in elem[1][5]:
					outfile.write("%s,%f\n"%(elem[0],k))

	outfile.close()
	return Order

def ListBasePairsFromStruct(Struct):  # return dic={structure:[liste de pairs de base ],....}
    lista = []
    stack = []
    ListUmp =[]
    for i in range(len(Struct)):  # sequence length
	if  Struct[i] ==".":
	    ListUmp.append(i)
        if Struct[i] == '(':  # Opening base pair
            stack.append(i)
        elif Struct[i] == ')':  # Closing base pair
            k = stack.pop()
            lista.append((k, i))#0-based
	
    return lista,ListUmp

def GetPseudoknots(Struct):
    lista = []
    stack = []
    for i in range(len(Struct)):  # sequence length
        if Struct[i] == '<':  # Opening base pair
            stack.append(i)
        elif Struct[i] == '>':  # Closing base pair
            k = stack.pop()
            lista.append((k, i))#0-based
    return lista

def classification(mutationrate):
	R="Unkown"
	if mutationrate<0.02:
		R="L"
	if 0.02<=mutationrate<=0.04:
		R="M"
	if mutationrate>0.04:
		R="H"
	return R

if __name__=="__main__":
	#python2.7 ComputeCorrelationsMaximalhelicesElimination.py --file Incell/Modified/5S_INCELL1_5SFREQUENCY_COUNTex.txt  --output 2020-10-21 --fasta alignment_sequences/5S.fa --amplicon 20 --p /home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2020-10-20/Profiling/5SProfiling.csv
	setofhelices=[]
	#EmbeddedHelices==False
	MinLenHelix=3
	File=args.file #-f
	OutputFolder=args.output #-o
	Fastafile=args.fasta #-R
	AmpliconF=args.amplicon #-a #20
	AmpliconB=args.backwardamplicon #-b #20
	Profile=args.profile #-p #maximal helix classes
	#RNAprofile=args.RNAprofile #-N # profiling output from  a sample
	nbrofheliceschosenhelices= args.Helixnbr #nbr of helices to be considered in cascading algorithm 
	#print Profile, RNAprofile
	if not os.path.exists(OutputFolder):
    		os.makedirs(OutputFolder)
    	RNA=File.split("/")[-1].split(".")[0]
	if 1==1:
		
                #outputIndexation=os.path.join(OutputFolder,RNA+"_MutationIndex.csv")
		#print outputIndexation
		######### Parse RNA sequence and structure
		fp=open(Fastafile)
		for i, line in enumerate(fp):
			if i==1:
				seq=line
			if i==2:
				SS=line
		fp.close()
                lenRNA=len(seq)
		#print lenRNA
		try:
   			SS
		except NameError:# #structure not provided
			SS="".join(["na" for i in range(lenRNA)])	
		#parse MM UU MU UM file		
		MM = np.zeros((lenRNA,lenRNA)) #MM dict
		UU= np.zeros((lenRNA,lenRNA)) # uU dict
                UM= np.zeros((lenRNA,lenRNA))
		MU= np.zeros((lenRNA,lenRNA))
		N= np.zeros((lenRNA,lenRNA)) 
		Prior=np.zeros((lenRNA,lenRNA)) 
		ListPairs=ListBasePairsFromStruct(SS)[0]+ GetPseudoknots(SS)
		ListUnpaired=ListBasePairsFromStruct(SS)[1]
				

		with open(File) as f:  #i-1based 	 j-1based 	 UU 	 UM 	 MU 	 MM 	 #Reads 	 MM/#Reads  	 UU/#Reads 
		    next(f)
		    for line in f:
		       #(i,j,bb,kk,dd,aa,X,Y,Z)= line.split() # format with 9 columns 
		       #(i,j,bb,kk,dd,aa)= line.split()
		       (i,j,bb,kk,dd,aa,X)= line.split() # format with 7 columns 
		       # change to 0-based
		       I=int(i)-1
		       J=int(j)-1
		       #print I,J
		       #print I,J,int(MM)
		       #print line
		       #print i,j,I,J
		       #print bb
		       N[I,J]=int(bb+kk+dd+aa)
		       
		       MM[I,J]= int(aa)
		       UU[I,J]= int(bb) #0-based
		       UM[I,J]=int(kk)
		       MU[I,J]= int(dd)
		'''
		with open (outputIndexation, 'w') as outfile:	
			outfile.write("%s,%s,%s,%s,%s\n"% ("position","nucleotide","structure","MutationRate","MutationRange"))
			
		'''
		if 1==1:				
			for i in range(AmpliconF, lenRNA-AmpliconB):
				if i in ListUnpaired:
					L="Up"
				else:
					if i+1 not in  ListUnpaired and i-1 not in ListUnpaired:
						L="P"
					else:
						L="HE"
				X=0
				if (UU[i,i]+MM[i,i])!=0:
					X=MM[i,i]/float(UU[i,i]+MM[i,i])
					
				#outfile.write("%i,%s,%s,%.4f,%s\n"% (i+1,seq[i],L,X,classification(X)))

				####################### Compute Prior on the range [amplicon,lenRNA-amplicon]x[i+7 , lenRNA-amplicon]: 
				for j in range(i+7, lenRNA-AmpliconB):
					Totalcount=UU[i,j]+MM[i,j]+UM[i,j]+MU[i,j]
					if (float(UU[i,i]+MM[i,i])!=0 and float(UU[j,j]+MM[j,j])!=0):
						Prior[i,j]=MM[i,i]/float(UU[i,i]+MM[i,i])*MM[j,j]/float(UU[j,j]+MM[j,j])*Totalcount
					

		
		######### Parse maximal helices file
		with open(Profile) as csv_file:
			    csv_reader = csv.reader(csv_file, delimiter=',')
			    line_count = 0
			    for row in csv_reader:
				if line_count ==0:
				    line_count += 1
				else:
				    setofhelices.append((row[0],(int(row[1]),int(row[2]),int(row[3]))))				    

		csv_file.close()

		Thermodynamic=[]
		####### add a third feature that is telling about the thermodynamic stability, step1 initialize it with "Not sampled"
		setofhelices=[(elem[0],elem[1],"NotSampled") for elem in setofhelices]
		###### Step2 assign values from thermodynamic profiling
		for k in range(len(setofhelices)):
			for match in Thermodynamic:
				#print setofhelices[k][1],match[0]
				if setofhelices[k][1]==match[0]:
					lst = list(setofhelices[k])
					lst[2]=match[1]
					setofhelices[k]=tuple(lst)					
		######## Compute Rd
		OrderedMH=Computecorrelations(MM,UU,Prior,ListPairs,lenRNA,File,MinLenHelix,setofhelices)
		# prepare helices for the prediction
		OrderedhelicesWithCoord=[]
		for helix in OrderedMH:
			for helix2 in setofhelices:
				if helix==helix2[0]:
					OrderedhelicesWithCoord.append( helix2 )
					
		#print OrderedhelicesWithCoord
		print ("Prior and Rd successfully computed for "+(File.split("/")[-1]).split("_")[0]+ " RNA for the condition "+ (File.split("/")[-1]).split("_")[1] ) 	
