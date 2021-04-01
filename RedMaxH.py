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
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
#from yellowbrick.cluster import KElbowVisualizer
from sklearn.datasets import make_blobs
regex = re.compile('[^a-zA-Z]')
#First parameter is the replacement, second parameter is your input string
#regex.sub('', 'ab3d*E')
#Out: 'abdE'
parser = argparse.ArgumentParser()                                               
parser.add_argument("--file", "-f", type=str, required=True)
parser.add_argument("--profile", "-p", type=str, required=True)
#parser.add_argument("--RNAprofile", "-N", type=str, required=True)
parser.add_argument("--output", "-o", type=str, required=True)
parser.add_argument("--fasta", "-R", type=str, required=True)
parser.add_argument("--amplicon", "-a",  type=int, required=True)
parser.add_argument("--backwardamplicon", "-b",  type=int, required=True)
parser.add_argument("--Helixnbr", "-x",  type=int, required=True)
args = parser.parse_args()


def Merge(dict1, dict2):
    return(dict2.update(dict1))
    
def Computecorrelations(MM,UU,Prior,Pairs,lenRNA,Name,MinLenHelix,setofhelices):
	Norm= np.zeros((lenRNA,lenRNA)) 
	for i in range( len(MM)):
		for j in range(len(MM)):
			if Prior[i,j]!=0:
				Norm[i,j]=(MM[i,j]-Prior[i,j])/float(Prior[i,j])
				print "mmmmmmmmmm",Norm[i,j]
	#set diagonals to Zero
        for i in range( len(MM)):
		MM[i,i]=np.nan
		Norm[i,i]=0
	
	Norm[Norm==0]=np.nan # set null normalized values to NA because Norm was initially defined as matrix of zeros
	######### Output MM count and Rltdiff
	with open (os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+"_Rawdata.csv"), 'w') as outfile: 
		outfile.write("%s,%s,%s,%s,%s\n"% ( "i","j","ncts","MM","Rltdiff" ))
		for i in range( AmpliconF, lenRNA-AmpliconB):
			for j in range( i+7,lenRNA-AmpliconB):
				#print "hhhhhhhhhhhhhhhh",lenRNA-AmpliconB
				outfile.write("%i,%i,%s,%.4f,%.4f \n"% (i+1,j+1,str(seq[i]+seq[j]),round(MM[i,j],4),round(Norm[i,j],4)))
	outfile.close()
	######## output data for specific helices
	outputcorrelation=os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+"_Norm_helices.csv")
	outputcorrelation1=os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+"_RedmaxH_Rankedhelices.csv")
	contributingpairs=os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+"_contributingpairsforrankedhelices.csv")
	Lista=[]
	print "pairs in helices with length L\n" 
	for Helix in setofhelices:
		if Helix[1][2]>MinLenHelix: 	
			listtuple=[]	
			i=Helix[1][0]-1
			j=Helix[1][1]-1
			l=Helix[1][2]
			for k in range(l):# neigborhood [i-1,i+1]x[j-1,j+1]
				for (n,t) in [(i+k,j-k),(i+k,j-k-1),(i+k,j-k+1),(i+k-1,j-k),(i+k-1,j-k-1),(i+k-1,j-k+1),(i+k+1,j-k),(i+k+1,j-k-1),(i+k+1,j-k+1)]:
					if j+1<lenRNA and math.isnan(Norm[n,t])!=True and (n,t,Norm[n,t])not in listtuple:  
						Lista.append((Helix[0],n,t,Norm[n,t],l))
						#Listadit.append(((Helix[0],l),n,t,Norm[n,t]))
						listtuple.append((n,t,Norm[n,t]))
				#for (n,t) in [(i+k,j-k)]:
				#	print n+1,t+1,l		
        # Get the #of common pairs between helices NA not included.
	dicto= defaultdict()
	dictobis= defaultdict() #oneadditional element that is the length of the helix
	CommonPairs=defaultdict()
	Matrixcommon=defaultdict(lambda: defaultdict(lambda:0))
	#SortedAverageHelix=defaultdict()
	#initialization
	for elem in Lista:
		dicto[elem[0]]=list()
		dictobis[elem[0]]=list()
		CommonPairs[elem[0]]=list()	
	#incrementation
	for elem in Lista:
		dicto[elem[0]].append((elem[1],elem[2],elem[3]))
		dictobis[elem[0]].append((elem[1],elem[2],elem[3],elem[4]))
	#print "ddddddddddddd",dicto
	Averagehelix=defaultdict() 
	Averagehelixtoadd=defaultdict() 
	dict3=defaultdict()
	for key in dictobis:
		
		toaverage=[]
		## CC content and AA content 
		listAA=0
		listCC=0
		for elem in dictobis[key]:
			# verification if key=="eP19plus1BP":
			#	print elem, seq[elem[0]],seq[elem[1]]
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
	
	
	'''
	#1 ----------first filtering layer not enough data, the nbr of observations is strctely less than the half of the expected number of observations 
	for elem in Averagehelix.keys():
		
		if Averagehelix[elem][1]<Averagehelix[elem][2]/2:
			del Averagehelix[elem]
	#print "After first filtering", Averagehelix
	#print Averagehelix
	'''
	#2 ------------ look for common pairs and remove those that are dominated by other helices
	for key1 in Averagehelix:
		firstlist=[elem for elem in dicto[key1]]
		for key2 in Averagehelix:
			if key2!=key1:
				secondlist=[elem for elem in dicto[key2]]
				#get pairs in common btw helices
				
				if len(set(firstlist) & set(secondlist))!=0:
					CommonPairs[key1].append((key2,len(set(firstlist) & set(secondlist))))
				#print "ddddddddddddd",key1,key2,len(set(firstlist) & set(secondlist))
				Matrixcommon[key1][key2]=len(set(firstlist) & set(secondlist))
				#print [elem for elem in dicto.keys()]
	#print  set(firstlist) ,"\n",  set(secondlist),"\n", set(firstlist) & set(secondlist)
	#print CommonPairs	
	removed=[]
	#print [ i for i in Matrixcommon.keys()]#:
	#for j in  Matrixcommon[i].keys():
	#		print i,j,Matrixcommon[i][j],"\t"
	#	print "\n"
	#print Matrixcommon['C1']['C1']
	Commonpairsfile="2DHelicesCommonPairs.csv"
	with open (os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+Commonpairsfile), 'w') as outfile:
		outfile.write("%s\n"% ( str([Helix[0] for Helix in setofhelices])))#,"commonpairs" ))
		for i in  [Helix[0] for Helix in setofhelices]:
			outfile.write("%s"% (i))
			for j in [Helix[0] for Helix in setofhelices]:
				#print i,j
				outfile.write(",%i"% (Matrixcommon[str(i)][str(j)]))
			outfile.write("\n")
		outfile.close()
	Commonpairsfile="MergedHelices.csv"
	with open (os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+Commonpairsfile), 'w') as outfile:	
		
	        outfile.write("%s\t%s\t%s\t%s\n"% ("Helix 1","Helix 2","# common pairs","Rc"))
		#combine instead of eliminate
		Donehelices=[]
		for Helix in Averagehelix:
			for Helixrelated in CommonPairs[ Helix]:
				X=[]
				#outfile.write("%s%s"% (Helix,Helixrelated))
				Helice=Helixrelated[0]	
				commonpairs=Helixrelated[1]
				
		
				# if dominated by another helix in average and most than half of observed values are in common, remove this helix
				#print Averagehelix[Helice][0],Averagehelix[Helix][0],commonpairs,Averagehelix[Helix][1],Averagehelix[Helice][0]<Averagehelix[Helix][0] , commonpairs>=Averagehelix[Helix][1]
				#get the list of all pairs
				#get pairs in common btw helices
				
				#print Helix, Helice, Averagehelix[Helix][5],Averagehelix[Helice][5],list(set([elem for elem in dicto[Helix]]) & set([elem for elem in dicto[Helice]]))
				#print "origi",Averagehelix[Helix][5], Averagehelix[Helice][5]
				X=Averagehelix[Helix][5]
				#print "fffffffffffff",Helix,Helice,len(Averagehelix[Helix][5]),len(Averagehelix[Helice][5]),len(set([elem for elem in dicto[Helix]]) & set([elem for elem in dicto[Helice]]))
				X= Averagehelix[Helix][5]+Averagehelix[Helice][5]
				#print X,"beforelll",list(set([elem for elem in dicto[Helix]]) & set([elem for elem in dicto[Helice]]))
				for elem in list(set([elem for elem in dicto[Helix]]) & set([elem for elem in dicto[Helice]])):
					if elem[2]>3:
						X.remove(3)
					else:
						X.remove(elem[2]) # remove will romove one instance of a value
				#print "afterafterafetr",len(X)
					
				if (Helix,Helice) not in Donehelices:
					Donehelices.append((Helix,Helice))
					Donehelices.append((Helice,Helix))
					if  commonpairs>9:#minimum 10 un BP moore neigborhoood:
						outfile.write("%s\t%s\t%i\t%.3f\n"% (Helix,Helice,commonpairs,round(np.mean(X),3)))
						Averagehelixtoadd[str(Helix+"-"+Helice)]=(round(np.mean(X),3),len(X),999,999,999,X,round(np.std(X),3),round(np.median(toaverage),3),round(np.mean([x for x in toaverage if x>0]),3),round(np.std([x for x in toaverage if x>0]),3),round(np.median([x for x in toaverage if x>0]),3))
	
						#feed the list of Rc with the new combinaisons and remove the preexisting ones.
						print Helix,Helice,X
						removed.append(Helix)
						removed.append(Helice)
				#if Averagehelix[Helix][0]<Averagehelix[Helice][0] and 2*commonpairs>=Averagehelix[Helix][1]:
				#print "removed ", Helix, "by whom", Helice
				#
	outfile.close()
	for Helix in list(set(removed)): # a trick to remove redundant helices
		#print Helix,"kkk",Averagehelix[Helix]
		if Averagehelix[Helix]:
			del Averagehelix[Helix]#
			#print "this is what we will remove",Averagehelix[Helix]
			# (-0.148, 21, 24, 0, 0, [-1.0, -1.0, -1.0, -0.5133488844147525, 0.12864787639752948, -0.05882939069427924,...
			
	#add combined helices 
	
	#remove helices with legth <4:
	shorthelices=[]
	
	for helix in Averagehelix:
		print key
		if Averagehelix[helix][1]<(3*5+4):
			shorthelices.append(helix)
	for helix in list(set(	shorthelices)):
		del Averagehelix[helix]
	#Averagehelix[key]=(round(np.mean(toaverage),3),len(toaverage),dictobis[key][0][3]*5+4,listCC,listAA,toaverage)
	#print "all",Averagehelix,Averagehelixtoadd,
	# combine two sets
	for k, v in chain(Averagehelix.items(),Averagehelixtoadd.items()):
		#print k,v
    		dict3[k]=v
 	#print dict3
	#sort the dictionary:
	SortedAverageHelix = dict3.items()
	SortedAverageHelix.sort(key=lambda x:-x[1][0])
	#print  "initila,",SortedAverageHelix
	#print "Averagehelix \n",Averagehelix	
			
	with open (outputcorrelation, 'w') as outfile:
		outfile.write("%s,%s,%s,%s,%s,%s,%s\n"% ( "Dataset","Helix","i","j","seqij","Rltdiff","ijdist"))#,"commonpairs" ))
		for elem in Lista:
			outfile.write("%s,%s,%s,%s,%s,%.4f,%i \n"% ( ((Name.split("/")[-1]).split("_")[1]).split(".")[0],elem[0],str(elem[1]+1),str(elem[2]+1),str(seq[elem[1]]+seq[elem[2]]),round(elem[3],4),elem[2]-elem[1]+1))#,CommonPairs[elem[0]]))
	
	Order=[]
	with open (outputcorrelation1, 'w') as outfile:
		#here the number of clusters
		outfile.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"% ( "Helix","averageRLtdiff","nbrobservations","Length","nbrCCpairs","nbrAApairs","Profiling","std","median","averageonpositive","stdonpositive","medianonpsitive" ))
		for elem in SortedAverageHelix:
			Order.append(elem[0])
			#get profiling
			if '-' not in elem[0]:
				outfile.write("%s,%.4f,%i,%i,%i,%i,%s,%.4f,%.4f,%.4f,%.4f,%.4f\n"% (elem[0], elem[1][0] , elem[1][1], (elem[1][2]-4)/5, elem[1][3], elem[1][4], [i[2] for i in setofhelices if i[0]==elem[0]][0], elem[1][6], elem[1][7], elem[1][8], elem[1][9], elem[1][10] ))
			else:
				# a combined profile
				outfile.write("%s,%.4f,%i,%i,%i,%i,%s,%.4f,%.4f,%.4f,%.4f,%.4f\n"% (elem[0], elem[1][0] , elem[1][1], (elem[1][2]-4)/5, elem[1][3], elem[1][4],  [i[2] for i in setofhelices if i[0]==elem[0].split("-")[0]][0] + " &  "+ [i[2] for i in setofhelices if i[0]==elem[0].split("-")[1]][0] , elem[1][6], elem[1][7], elem[1][8], elem[1][9], elem[1][10]))
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
	
	("P3a",(21,68,8))
	
def EmbeddedHelices(listofHelices): # Embedded helices are defined as follow: H2 is inside the space [H1+L1,H1-L1] the maximal distance is two nucleotides from each side.
	List=[]
	#listofHelices=[("Cx",(79,97,8)),("Cx2",(83,92,4)),("Cx3",(80,98,4))]
	for H1 in listofHelices:
		for H2 in listofHelices:
			#print  H1,H2,H2[1][0]>H1[1][0],H2[1][0]-(H1[1][0]+H1[1][2])<3,H1[1][1]-H1[1][2]-H2[1][1]<3
			if H2[1][0]>=H1[1][0]+ H1[1][2] and H2[1][1]<=H1[1][1]- H1[1][2]and  H2[1][0]-H1[1][0]-H1[1][2]<3 and H1[1][1]-H1[1][2]-H2[1][1]<3:
				# length of the new embeededhelix is the sum of the two helices lengths and the maximal gap value that is less than 3
				LengthEmbeddedhelix=H1[1][2]+H2[1][2]+np.amax([H2[1][0]-H1[1][0]-H1[1][2],H1[1][1]-H1[1][2]-H2[1][1]])
				List.append((str(H1[0]+H2[0]),(H1[1][0],H1[1][1],LengthEmbeddedhelix)))
	#print "List embeeded", List
	return List

def RunEval(InputFile):
    Energy=0
    # launch the RNaeval command
    #conf = loadConfig()
    energiesFile =  "energyvalues.txt"
    os.system('RNAeval -d2 <' + InputFile  + '>' + energiesFile) # -v if we need all loops energy values
    # Parse the RNAevaloutput to extract energy values
    lines = PR.Parsefile(energiesFile)
    #print lines
    for i in xrange(1, len(lines), 2):
        # i is the stucture number and 'lines[i].split(" ")[1][1:-2]' is  the  corresponding  energy value
        Energy= (lines[i].split(" ")[1][1:-2])
    return Energy
       
if __name__=="__main__":
	#python2.7 ComputeCorrelationsMaximalhelicesElimination.py --file Incell/Modified/5S_INCELL1_5SFREQUENCY_COUNTex.txt  --output 2020-10-21 --fasta alignment_sequences/5S.fa --amplicon 20 --p /home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2020-10-20/Profiling/5SProfiling.csv
	
	

	setofhelices=[]
	EmbeddedHelices==False
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
                lenRNA=len(seq)-1
		
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
		       print I,J
		       #print I,J,int(MM)
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
					
		#outfile.close()
		
		######### Parse maximal helices file
		with open(Profile) as csv_file:
			    csv_reader = csv.reader(csv_file, delimiter=',')
			    line_count = 0
			    for row in csv_reader:
				if line_count ==0:
				    line_count += 1
				else:
				    setofhelices.append((row[0],(int(row[1]),int(row[2]),int(row[3]))))				    
		#print "onee,",setofhelices
		csv_file.close()
		#print "twoo", EmbeddedHelices(setofhelices)
		#if EmbeddedHelices:
		#setofhelices=setofhelices+EmbeddedHelices(setofhelices)
		#print setofhelices
		Thermodynamic=[]
		######### Parse Profiles generate from a boltzmann sample
		#with open(RNAprofile) as csv_file:
		#	    csv_reader = csv.reader(csv_file, delimiter=',')
		#	    line_count = 0
		#	    for row in csv_reader:
		#		if line_count ==0:
		#		    line_count += 1
		#		else:
				   
				   	 #Thermodynamic.append(((int(row[0]),int(row[1]),int(row[2])),row[3]))		
		#csv_file.close()
		####### add a third feature that is telling about the thermodynamic stability, step1 initialize it with "Not sampled"
		setofhelices=[(elem[0],elem[1],"NotSampled") for elem in setofhelices]
		###### Step2 assign values from thermodynamic profiling
		for k in range(len(setofhelices)):
			for match in Thermodynamic:
				#print setofhelices[k][1],match[0]
				if setofhelices[k][1]==match[0]:
					#print here ,setofhelices[k][2] 
					#print setofhelices[k][2],match[1] a conversion to a list seems to be necessary in order to update the value Notsampled
					lst = list(setofhelices[k])
					lst[2]=match[1]
					setofhelices[k]=tuple(lst)
					#print "after" , setofhelices[k][2]
					
		##### for all helices generate the profile
		#for elem in setofhelices:
		#	print "all",elem
		######## Compute RLtdiff
		OrderedMH=Computecorrelations(MM,UU,Prior,ListPairs,lenRNA,File,MinLenHelix,setofhelices)
		# prepare helices for the prediction
		OrderedhelicesWithCoord=[]
		for helix in OrderedMH:
			for helix2 in setofhelices:
				if helix==helix2[0]:
					OrderedhelicesWithCoord.append( helix2 )
					
		#print OrderedhelicesWithCoord
		print ("Prior and Rltdiff successfully computed for "+(File.split("/")[-1]).split("_")[0]+ " RNA for the condition "+ (File.split("/")[-1]).split("_")[1] ) 
	'''
	### Structure prediction using the set of sorted helices based on a decreasing Rltdiff average value
			
	seq=seq.strip() # clean the sequence from \n
	SS=SS.strip() #idem pour la structure
	print "Native structure",SS
	Method=4#4 1 2 or 3 5
	#nbrofheliceschosenhelices=11 # how many HC is going to be tested
	##combi=30 # number of combination of constraints
	Predictedstructures=[]
	PredictedstructuresCombinations=[]
	#listofhelices=OrderedhelicesWithCoord
	# This transformatin is necessary to keep using tuple but in a list
	listofhelices1=[[elem] for elem in OrderedhelicesWithCoord]
	Nativehelices=[[elem for elem in setofhelices if elem[0][0]!="C"]] #all helices
	#print "which helix was ssleected" ,OrderedhelicesWithCoord
	#print "hhhhh", Nativehelices
	#Combin=PR.combinations(listofhelices[:nbrofheliceschosenhelices])
	#print Combin ,"before"
	#Combin=PR.ComputeCombinations(listofhelices[:nbrofheliceschosenhelices])
	#print listofhelices1
	Combin=PR.CascadingConstraintsversion2(OrderedhelicesWithCoord[:nbrofheliceschosenhelices])
	#print Combin ,"after"fv
	Combin=[[elem] for elem in Combin]
	#print "hhh", listofhelices[:nbrofheliceschosenhelices]
	#########1-MFE
	Predictedstructures.append(PR.Predict2D(Fastafile,"MFE"))
	#######2- Elements in the sorted helices set 
	for case  in listofhelices1:
		#print "case",case
		Constraintfile=PR.buildConstraint(seq, case,Method)
		#print "fff",Constraintfile
		#print  PR.Predict2D(Constraintfile,str(case)+str(Method))
		Predictedstructures.append( PR.Predict2D(Constraintfile,str(case)+str(Method)))
	########3- Combinations 
	
	#Combin[:0] =Nativehelices # append native helices at the beginning, pay attention because there is an overlaping issue!
	Combin[:0] =PR.CascadingConstraints(Nativehelices)
	#print "hhhhh"
	#print "eeeeeeeeeee",Nativehelices,"Possible combinaions =" ,Combin
	#print "llllllllllll",Combin
	combi=len(Combin)
	
	for X in range(combi):
		for case in Combin[X]:#listofhelices1+Combin[1:]:
			#print case
			#print "hopla",case
			Constraintfile=PR.buildConstraint(seq, case,Method)#
			PredictedstructuresCombinations.append( PR.Predict2D(Constraintfile,str(case)+str(Method)))
		
	
	##############Cleaning all ps RNAfold generated plots 
	dir="."	
	PR.removefileswithextension(dir)
	##########Create a structure folder
        StructureFolder="./"+ str(RNA)
        if not os.path.exists(StructureFolder):
    		os.makedirs(StructureFolder)
	#################Assess the prediction accuraccy
	BP=  PR.ListBasePairsFromStruct(SS)
	#sprint "pairs from the refernce", BP
	outputaccuracy=os.path.join(OutputFolder,(RNA+"_accuracyMaximalHelixConstrainedStructure.csv"))
	outputaccuracyCombinations=os.path.join(OutputFolder,(RNA+"_accuracyMaximalHelixConstrainedStructureCombinations.csv"))
	Ind=0
	with open (outputaccuracy, 'w') as outfile:
		outfile.write("%s,%s,%s,%s,%s\n"% ( "Constraint","Sensitivity","PPV","Accuracy","Structure" ))
		for PredictedStr in Predictedstructures:
			BPM= PR.ListBasePairsFromStruct(PredictedStr[2].split(" ")[0]	)
			MCCSHAPE_Rsample=(PR.PPV_sensitivity_Rsample(BP,BPM)[0],PR.PPV_sensitivity_Rsample(BP,BPM)[1])
			outfile.write("%s,%.4f,%.4f,%.4f ,%s\n"% (PredictedStr[0].replace(',', ''), round(MCCSHAPE_Rsample[0],4), round(MCCSHAPE_Rsample[1],4), round(math.sqrt(round(MCCSHAPE_Rsample[0],4)*round(MCCSHAPE_Rsample[1],4)),3),PredictedStr[2].split(" ")[0]))
			if Ind==0:
				###"get the accuracy of the MFE
				rapport0=math.sqrt(round(MCCSHAPE_Rsample[0],4)*round(MCCSHAPE_Rsample[1],4))
				print rapport0
			Ind+=1
	outfile.close()
	
	Pstructures=[]
	with open (outputaccuracyCombinations, 'w') as outfileX:
		outfileX.write("%s,%s,%s,%s,%s,%s,%s,%s\n"% ( "Constraint","Sensitivity","PPV","Accuracy","%Accuracy/optimal","%Accuracy/MFE","energy","Structure" ))
		Index=0
		#print PredictedstructuresCombinations
		for PredictedStr in PredictedstructuresCombinations:
			
			#### Generate structures to visualize
			#os.makedirs(OutputFolder+"/"+RNA)
			outfilestructure=  os.path.join(StructureFolder,(PredictedStr[0]+ ".dbn"))
			with open (outfilestructure  , 'w') as outfilest:
				outfilest.write("%s\n%s\n%s"% (PredictedStr[0].replace(',', ''),seq,PredictedStr[2].split(" ")[0]))
			outfilest.close()
			
			##########Generate Files for RNAeval ##keep this order dbn than txt
			outfilestructure=  os.path.join(StructureFolder,(regex.sub('', PredictedStr[0])+ ".txt"))
			with open (outfilestructure  , 'w') as outfilest:
				outfilest.write("%s\n%s"% (seq,PredictedStr[2].split(" ")[0]))
			outfilest.close()
			######################compute energy 
						
			BPM= PR.ListBasePairsFromStruct(PredictedStr[2].split(" ")[0]	)
			#print PredictedStr[2]
			MCCSHAPE_Rsample=(PR.PPV_sensitivity_Rsample(BP,BPM)[0],PR.PPV_sensitivity_Rsample(BP,BPM)[1])
			
			if Index==0:###"get the accuracy of the reference structure
				rapport=math.sqrt(round(MCCSHAPE_Rsample[0],4)*round(MCCSHAPE_Rsample[1],4))
				#print rapport
			outfileX.write("%s,%.4f,%.4f,%.4f,%.2f,%.2f,%.1f,%s\n"% (PredictedStr[0].replace(',', ''), round(MCCSHAPE_Rsample[0],4), round(MCCSHAPE_Rsample[1],4), round(math.sqrt(round(MCCSHAPE_Rsample[0],4)*round(MCCSHAPE_Rsample[1],4)),3),math.sqrt(MCCSHAPE_Rsample[0]*MCCSHAPE_Rsample[1])/float(rapport)*100,math.sqrt(MCCSHAPE_Rsample[0]*MCCSHAPE_Rsample[1])/float(rapport0)*100,float(RunEval(outfilestructure)),PredictedStr[2].split(" ")[0]))
			
			Pstructures.append(PredictedStr[2].split(" ")[0])
			Index+=1	
			#############
		outfileX.write("%s %s \n"% ("# different structures = ", len(set(Pstructures))))		
	outfileX.close()
	'''		
