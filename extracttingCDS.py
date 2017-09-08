import re
### reading in file using the readTestFile function
def readTestFile(file):
	handle=open(file)
	lines="\n"
	for line in handle:
		lines+=line
	handle.close()
	return lines
	
def getORF(seqToRead):
	orfs=[]
	orfs=re.findall('((ATG){1}(\w{3})+(TAG|TGA|TAA){1})',seqToRead)
	return orfs

def getonlyFirstElement(tupleOfSequences):
	onlyORFs=[]
	for listOfSeq in tupleOfSequences:
		orf=listOfSeq[0]
		##print(orf)
		onlyORFs.append(orf)
	return onlyORFs
	
def getORFWithHighestLength(listOfORFs):
	ORFWithHighestLen=""
	ORFWithHighestLen=max(listOfORFs,key=len)
	return ORFWithHighestLen

file='BRCA2.fa'
##print(readTestFile(file))
fastaRead=readTestFile(file)
listOfSeq =getORF(fastaRead)
allORFsRetrieved = getonlyFirstElement(listOfSeq)
##print(allORFsRetrieved)
print("The ORFs with highest length is")
print(getORFWithHighestLength(allORFsRetrieved))
print(len(getORFWithHighestLength(allORFsRetrieved)))		
 