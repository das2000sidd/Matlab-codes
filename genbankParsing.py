import re
## reading in file using readtTestFile
def readTestFile(file):
	handle=open(file)
	lines="\n"
	for line in handle:
		lines+=line
	handle.close()
	return lines


### retrieving get sections in file lines read
def getMainSections(file1):
	
	segments = re.split("\n(\w+)",file1)
	map={}
	for i in range(1,len(segments),2):
		segments[i]=segments[i].rstrip('\r\n')
		map[segments[i]]=segments[i+1]
	return map

string="spo11.gb"
fileRead=readTestFile(string)
testingDict = getMainSections(fileRead)
print("The dictionary of all features")
##print(testingDict) ## WORKING FINE
##print(testingDict) 


## retrieving all protein IDs
def findAllProteinID(file2):
	string=readTestFile(file2)
	proteinIDs = re.findall("([X|N]P_\d+\.\d)",string) ### protein ID start with X or N, followed by a "P_" and then several numbers, then a dot and then one number
	return proteinIDs

allProteinIDs=findAllProteinID(string)
print(allProteinIDs) ### working fine


### get start of exons from sequence titled CDS

def getAllExonLoci(dictionaryOfFeatures):
	allCDSLoci=[]
	for key,value in dictionaryOfFeatures.items():
		if (key=='FEATURES'):
			 
			value=value.replace('"','') ## getting rid of quotes
			##print(value)
			## logic of finding the regex was considering each position at a time and then quantifying
			## for CDS start with many spaces and hence \s+
			## Then the word CDS comes and new regex is \s+(CDS)
			## Then again comes a lot of spaces and so regex changes to \s+(CDS)\s+
			## After (spaces)CDS(spaces) comes the word join and so regex changes to \s+(CDS)\s+(join)
			## After join comes one bracket with no repitition. So regex changes to \s+(CDS)\s+(join)\(
			## Now repitation states and for all positions we can have either numbers or . or , or \s
			## and so regex changes to \s+(CDS)\s+(join)\(([\d.,\s]+)
			## After that is over we just have one position that is the closing bracket and so regex changes to
			## \s+(CDS)\s+(join)\(([\d.,\s]+)\)
			allCDSLoci = re.findall("\s+(CDS)\s+(join)\(([\d\.,\s]+)\)",value)
			

	return allCDSLoci
	
geneCDSLoci = getAllExonLoci(testingDict)
print("The CDS loci are")
print(geneCDSLoci)	## working fine

### get all the mRNA annotation
def getAllmRNAAnnotation(dictionaryOfFeatures):
	allmRNALoci=[]
	for key,value in dictionaryOfFeatures.items():
		if (key=='FEATURES'):
			value=value.replace('"','')
			allmRNALoci = re.findall("\s+(mRNA)\s+(join)\(([\d\.,\s]+)\)",value)
	return allmRNALoci

genemRNAloci=getAllmRNAAnnotation(testingDict)
print("The mRNA loci are")
print(genemRNAloci)


### Retrieve one of the mRNA annotation from the list of mRNA annotations and getting its start and stop positions
anyOnemRNA=genemRNAloci[1]
lociListOfAnyOnemRNA=anyOnemRNA[2]
anyOnemRNAStartPositions = re.findall(r'(\d+\.\.)+',lociListOfAnyOnemRNA)	
mRNAStartPositions=[]
for mRNAPos in anyOnemRNAStartPositions:
	
	anyOnemRNAStartPositionsOnlyNo=re.findall(r'(\d+)+',mRNAPos)
	mRNAStartPositions.extend(anyOnemRNAStartPositionsOnlyNo)

anyOnemRNAStopPositions=re.findall(r'(\.\.\d+)+',lociListOfAnyOnemRNA)
mRNAStopPositions=[]
for mRNAPos in anyOnemRNAStopPositions:
	anyOnemRNAStopPositionsOnlyNo = re.findall(r'(\d+)',mRNAPos)
	mRNAStopPositions.extend(anyOnemRNAStopPositionsOnlyNo)
	
## Extract the loci start and stop separately- Only start done but working
def geneCDSLociGetOnlyStart(geneCDSLoci):
	
	eachCDS=()
	cdsStartingPosDict={}
	eachCDSLocus=""
	eachCDSStartingPos=[]
	for i in range(0,len(geneCDSLoci),1):
		eachCDS=geneCDSLoci[i]
		eachCDSLocus=eachCDS[2] ### the starting..ending of exons
		eachCDSLocusSplit=re.findall(r'(\d+\.\.)+',eachCDSLocus)
		## re.findall helped pulled out the pattern "numbers followed by two dots repeated multiple times"
		
		print(eachCDSLocusSplit)
		
print("The CDS starting positions are")	
print(geneCDSLociGetOnlyStart(geneCDSLoci))

	
## Retrieving the dna sequence
def getTheDNASequence(dictionaryOfFeatures):
	dnaSequence=""
	sequenceToReturn=""
	for key,value in dictionaryOfFeatures.items():
		if (key=='ORIGIN'):
			value=value.replace('\n','')
			### This regex tells find a pattern with only a or t or g or c or space such that thet occur multiple times next to each other
			dnaSequence= re.findall("([atgc\s]+)",value) ## list like this [' ccgctcagaa agcgcgggaa aggcacgcag ccacgcccca agggcgcagc ctaggacagg       ', ' ggcttctgga gcttctggca gccgtctgcc ctcatggcct ttgcacctat ggggcccgag      ',]
			for eachSequence in dnaSequence:
				eachSequence=eachSequence.replace(" ","")
				sequenceToReturn=sequenceToReturn+eachSequence
	return sequenceToReturn


dnaSeqRetrieved=getTheDNASequence(testingDict)


def convertLociAsStringToNumeric(lociList):
	lociListAsInt=[]
	for anyItem in lociList:
		intVal=int(anyItem)
		lociListAsInt.append(intVal)
	return lociListAsInt
	
mRNAStartPosAsInt = convertLociAsStringToNumeric(mRNAStartPositions)
mRNAEndPosAsInt = convertLociAsStringToNumeric(mRNAStopPositions)
print("This list is start positions of any one mRNA")
print(mRNAStartPosAsInt) ## no quotes
print("This list is end positions of any one mRNA")
print(mRNAEndPosAsInt) ## no quotes

def getFourBasesAroundStartOfExon(dnaSequence,mRNAStartPositions):
	dictToStoreSequenceAroundStart={}
	for startPos in mRNAStartPositions:
		if(startPos==1):
			basesBeforeStart=""
		else:
			basesBeforeStart=dnaSequence[startPos-4:startPos] ### slicing
		basesAfterStart=dnaSequence[startPos+1:startPos+5]
		basesAroundStart=basesBeforeStart+"+"+basesAfterStart
		dictToStoreSequenceAroundStart[startPos]=basesAroundStart
	return dictToStoreSequenceAroundStart
			
			
def getFourBasesAroundEndOfExon(dnaSequence,mRNAStopPositions):
	dictToStoreSequenceAroundEndOfExon={}
	for stopPos in mRNAStopPositions:
		if(stopPos==len(dnaSequence)):
			basesAfterEnd=""
		else:
			basesAfterEnd=dnaSequence[stopPos+1:stopPos+5]
		basesBeforeEnd=dnaSequence[stopPos-4:stopPos]
		basesAroundEnd=basesBeforeEnd+"+"+basesAfterEnd
		dictToStoreSequenceAroundEndOfExon[stopPos]=basesAroundEnd
	return dictToStoreSequenceAroundEndOfExon

listOfSeqAroundStartOfExon=getFourBasesAroundStartOfExon(dnaSeqRetrieved,mRNAStartPosAsInt)
listOfSeqAroundEndOfExon=getFourBasesAroundEndOfExon(dnaSeqRetrieved,mRNAEndPosAsInt)
##print(dnaSeqRetrieved)
print("These are the sequences around start of exon")
print(listOfSeqAroundStartOfExon)
print("These are the sequences around end of exon")
print(listOfSeqAroundEndOfExon)
print("The dna sequence is")
print(dnaSeqRetrieved)
		
		
### ALL PRACTICAL PARTS ANSWERED###

	