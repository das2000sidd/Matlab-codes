import re
### reading in file using the readTestFile function
def readTestFile(file):
	handle=open(file)
	lines="\n"
	for line in handle:
		lines+=line
	handle.close()
	return lines
## creating a dictionary where seq id is key and actual sequence is value
def getSections(file1):
	string=readTestFile(file1)
	segments=re.split("(>\w+\d+\|\w+\d+)",string) ## splits the fasta into list of fasta id and fasta sequence. (>\w+\d+\|\w+\d+) precedence group for (>FBgn0021761|FBtr0340224) format
	map={}
	for i in range(1,len(segments),2):
		segments[i]=segments[i].rstrip('\r\n')
		map[segments[i]]=segments[i+1]
	return map


### creating function for retrieving arginine codon percentage for all codons
def getArgininePercent(dictionaryOfmRNASeq): ## must pass dictionary
	dictOfCodonCount={}
	for key,value in dictionaryOfmRNASeq.items():
			noOfArg=0
			value=value.replace('\n','')
			for i in range(0,len(value),3):
				aCodon=value[i:i+3]
				if(aCodon in ["AGA","AGG"]):
					noOfArg=noOfArg+1
					##dictOfCodonCount[key]=noOfArg/(len(value)/3) ## Made an indentation error, INDENTATION VERY IMPORTANT
			dictOfCodonCount[key]=noOfArg/(len(value)/3) ## Made an indentation error, INDENTATION VERY IMPORTANT		
	return dictOfCodonCount
	


			
			
chr2LproteinSeq="2L_CDS.txt"
codingSeqDict = getSections(chr2LproteinSeq) ## dictionary is returned		
dictOfArgPercent = getArgininePercent(codingSeqDict)
text_file1=open("arginine percentage of proteins.txt","w")
for key,value in dictOfArgPercent.items():
	text_file1.write("The arginine percentage of " + key + " is " + str(value))	
	text_file1.write('\n')
text_file1.close()
	
	
	
### WORKING FINE



### looking at the file and inputting those list for GO annotation there were some (5) with DNA binding###
###FBgn0010287,FBgn0032517,FBgn0086445,FBgn0003732,FBgn0003607

### trying the read in CDS, translate and then find arginine percentage
def convertCodontoDictionary(codonFile):
	codonDictionary={}
	string=readTestFile(codonFile)
	segments=re.split("\n",string)
	##print(segments)
	if('' in segments):
		segments.remove('')
	for segment in segments:
		codon,aminoAcid=re.split("\t",segment)
		codonDictionary[codon]=aminoAcid
	return codonDictionary
	##segments=re.split("\n",lines)
	##codon,aminoAcidCode=re.split('\s+',segments)
	##codon=codon.strip()
	##aminoAcidCode=aminoAcidCode.strip()
	##codonDictionary[codon]=aminoAcidCode
	##return codonDictionary
	
codonTableFile="codon_coding_table.txt"
codonDict = convertCodontoDictionary(codonTableFile) ## dictionary of codon and amino acid code
##print(codonDict)
def getListOfAllCodons(codonFile):
	codonList=[]
	string=readTestFile(codonFile)
	segments=re.split("\n",string)
	if('' in segments):
		segments.remove('')
	for segment in segments:
		codon,aminoAcid=re.split("\t",segment)
		codonList.append(codon)	
	return codonList
allCodonList=getListOfAllCodons(codonTableFile)	## List of all codons
##print("These are the CDS regions")
##print(cdsRegionsChr2) ## dictionary of seqid and sequence
#### TRANSLATION###
def makeProteinFromCodon(CdsDict,CodonDict):
	translatedProteinSeqDict={}
	translatedProteinSequence=""
	for Seqid,CodingSequence in CdsDict.items():
		CodingSequence=CodingSequence.replace('\n','')
		for i in range(0,len(CodingSequence),3):
			aCodon=CodingSequence[i:i+3]
			aCodon=aCodon.lower()
			for codon,aminoAcid in CodonDict.items():
				if(aCodon==codon):
					translatedProteinSequence+=str(aminoAcid)
		translatedProteinSeqDict[Seqid]=translatedProteinSequence
	return translatedProteinSeqDict

translatedProteinSequence=makeProteinFromCodon(codingSeqDict,codonDict)
##print("These are the translated amino acids corresponding to the CDS regions")
##print(translatedProteinSequence) ### working fine


def getGenesMatchingMotif(dictionaryOfmRNASeq):
	listOfGenesWithPattern=[]
	for key,value in dictionaryOfmRNASeq.items():
		value=value.replace('\n','')
		patternfound=re.findall("([QRL].[EDK][VAIL][AGK]..[LAVM]G[VIL][ST]..[TQA][VIL][SR][RK])",value,re.M|re.I)
		##patternfound=re.findall("(HHEE)",value,re.M|re.I)- DO NOT DELETE THIS LINE. THIS IS LOGIC FOR FINDING SEQUENCES THAT HAVE THE PATTERN
		if(len(patternfound) > 0):
			
			listOfGenesWithPattern.append(key)
	return listOfGenesWithPattern
	
geneWithThePattern=getGenesMatchingMotif(translatedProteinSequence)

test_file=open("gene_ids_with_pattern.txt","w") ### list of genes matching motif###

for string in geneWithThePattern:
	##geneId,transcriptId=string.split('|')
	##symbol,actualGene = geneId.split('>')
	test_file.write(string)
	test_file.write('\n')
	
test_file.close()
### Get percentage of basic and amino acids###
def getPercentageofAcidicandBasicAminoAcids(translatedProteinSequenceDict):
	
	
	acidicAminoAcids=[]
	basicAminoAcids=[]
	bothAminoAcidPercentage=[]
	aminoAcidPercentage={}
	for seqID,ProteinSeq in translatedProteinSequenceDict.items():
		ProteinSeq=ProteinSeq.replace('\n','')
		acidicAminoAcids = re.findall("[D|E|K]",ProteinSeq)
		basicAminoAcids=re.findall("[R|H]",ProteinSeq)
		acidicAminoAcidPercentage = (len(acidicAminoAcids)/len(ProteinSeq))*100
		basicAminoAcidPercentage = (len(basicAminoAcids)/len(ProteinSeq))*100
		bothAminoAcidPercentage.append(acidicAminoAcidPercentage)
		bothAminoAcidPercentage.append(basicAminoAcidPercentage)
		aminoAcidPercentage[seqID]=bothAminoAcidPercentage
		bothAminoAcidPercentage=[]
		
	return 	aminoAcidPercentage

aminoAcidPercentages=getPercentageofAcidicandBasicAminoAcids(translatedProteinSequence)
##print("In list below, first percentage is acidic and second is basic")
##print(aminoAcidPercentages)	
			
### Get all the codons for all coding sequences- WORKING FINE
def getListOfCodonForEachSequence(DNASeq):
	codonPercentageDictionary={}
	codonList=[]
	codonCount=0
	for seqID,CodingSequence in DNASeq.items():
		CodingSequence=CodingSequence.replace('\n','')
		for i in range(0,len(CodingSequence),3):
			aCodon=CodingSequence[i:i+3]
			aCodon=aCodon.lower()
			codonList.append(aCodon)
		codonPercentageDictionary[seqID]=codonList
	return codonPercentageDictionary

seqIdByCodon = getListOfCodonForEachSequence(codingSeqDict)
##print(seqIdByCodon)
## I have list of all codons
### get percentage of all codons by CDS
###EXTREMELY IMPORTANT CODE TO GET CODON PERCENT FOR ALL SEQUENCES BY SEQUENCE ID AND CODON
def getPercentageOfEachCodonPerCodingSequence(listOfCodonBySeqID,allPossibleCodons):
	codonCountPerCodon={}
	codonPercentPerSeq={}
	for seqId,codonList in listOfCodonBySeqID.items():
		for codon in allPossibleCodons:
			codonCount = codonList.count(codon)
			codonCountPerCodon[codon]=codonCount/len(codonList)
		codonPercentPerSeq[seqId]=codonCountPerCodon
	return codonPercentPerSeq
		
testingPerCentageOfCodingSeqPerCodon=getPercentageOfEachCodonPerCodingSequence(seqIdByCodon,allCodonList)		
##print(testingPerCentageOfCodingSeqPerCodon)					
					
##allCodonList=getListOfAllCodons(codonTableFile)				
			
		