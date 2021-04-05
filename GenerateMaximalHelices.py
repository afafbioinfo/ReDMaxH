
import os, subprocess,argparse

def Parsefile(Path):
    fileIn = open(Path, "r")
    lines = [l.strip() for l in fileIn.readlines()]
    fileIn.close()
    return lines


def parseArguments():

    parser = argparse.ArgumentParser(description = "Generate maximal helices")
    parser.add_argument('--Fasta', help='Path to fasta sequence file')
    parser.add_argument('--RNA', help='RNA name')
    parser.add_argument('--ForwardAmpliconSize', type=int, default=0, help='The size of the amplicon in the forward sense (default = 0)')
    parser.add_argument('--BackwardAmpliconSize', type=int, default=0, help='The size of the amplicon in the backward sense (default = 0)')
    
    parser.add_argument('--minHelixsize', type=int, default=4, help="""Minimum required helix size(default = 4). """)
    
    args = parser.parse_args()   
        
    return args
       
def GenerateMaximalhelices(seq, Hlength):

  	
	output="result.txt"
	Command = 'ExternalScripts/HelixEnumeration-master/HelixEnumeration-master/Build/src/HelixEnumeration -s '
	Command += seq 
  	Command += 'k' + str(Hlength) #it seems that the argument is not working
	subprocess.call(Command, stdout=open(output, 'wb'), stderr=open(os.devnull, 'w'), shell=True)
     	filesize = os.path.getsize(output)
	if filesize != 0:
		Temp=(Parsefile(output))
		os.remove(output)
		return Temp
    	else:
    		return "No generated helices"
	

if __name__ == '__main__':
	args = parseArguments()
	#verbal = True # afaf verbal to false
	verbal = False
	Tag=args.RNA
	ampliconforward=args.ForwardAmpliconSize
	ampliconbackward=args.BackwardAmpliconSize
	HelixKnownannotation=[]
	Hlength=args.minHelixsize 
	#print "ll",args
	fp=open(args.Fasta)
	for i, line in enumerate(fp):
		if i==1:
			seq=line
		if i==2:
			SS=line
	fp.close()

	helices=GenerateMaximalhelices(seq[ampliconforward:len(seq)-ampliconbackward],Hlength)
	Helix=[]
	# increment by one 
	for i in range(len(helices)):
		if int(helices[i].split(" ")[2])>=Hlength:
			Helix.append(("C"+str(len(Helix)+1),(int(helices[i].split(" ")[0])+1+ampliconforward,int(helices[i].split(" ")[1])+1+ampliconforward,int(helices[i].split(" ")[2]))))
			
	for elemx in range(len(Helix)):
		for elemy in HelixKnownannotation:
			if Helix[elemx][1]==elemy[1]:
				Helix[elemx]=elemy

	outputfile=os.path.join("RedMaxHoutput",Tag+"Profiling.csv")
	o = open(outputfile, "w")  
	o.write("%s,%s,%s,%s\n" % ("Helix","opening","closing","length"))
	for elem in Helix:
	    	o.write("%s,%i,%i,%i\n" % (elem[0], elem[1][0], elem[1][1],elem[1][2]))
	o.close()
