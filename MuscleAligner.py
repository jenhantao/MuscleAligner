#! /usr/bin/python

# Process: use muscle to do pairwise alignment of the reference to all of the short reads/sequences
# output: combine all of the aliignments and show coverage of the reference.  Then we can compute things like, where is there a gap in the coverage, and where do two different reads disagree?

# input file: fasta file, with first sequene being the reference
# second argument is whether or not to silence muscle, silenced by default, n to output everything
import os
import sys
import re

# reverse complements input strings
def revComp(s):
    s = s.upper()
    toReturn = ""
    for i in range (len(s)):
        if s[i] == "A":
            toReturn = "T" + toReturn
        elif s[i] == "T":
            toReturn = "A" + toReturn    
        elif s[i] == "C":
            toReturn = "G" + toReturn 
        elif s[i] == "G":
            toReturn = "C" + toReturn 
        else:
             toReturn = s[i] + toReturn
    return toReturn

if len(sys.argv) > 2:
    if sys.argv[2] == "n":
        quiet = ""
    else:
        quiet = "-quiet "
else:
    quiet = "-quiet "
with open("./"+sys.argv[1]) as f:
    inputSequenceFile = f.readlines()
sequenceNames = [];
sequences = [];
i=0
maxNameLength = 0 ;
while i<len(inputSequenceFile):
    line = inputSequenceFile[i]
    if line.startswith(">"):
        name = line[1:].strip()
        if len(name) > maxNameLength:
            maxNameLength = len(name)
        sequenceNames.append(name)
    else:
        sequence = ""
        while i < len(inputSequenceFile) and not line.startswith(">"):
            line = inputSequenceFile[i]
            if (not line.startswith(">")):
                sequence = sequence + line.strip() 
                i = i+1
            else:
                i = i-1
        sequences.append(sequence.strip())
    i = i + 1

# call muscle to do pairwise alignments
referenceSequence = sequences[0]
referenceName = sequenceNames[0]
i = 1
tempfile = open('temp.fas','w')
reads = []; # array that stores the alignments of the reads
references = []; # array that stores the alignments of the reference sequence
while i < len(sequenceNames):
    tempfile = open('temp.fas','w+')
    tempfile.write(">" + referenceName + "\n" + referenceSequence + "\n>" + sequenceNames[i] + "\n" + sequences[i] + "\n>" +sequenceNames[i] + "(-)\n" + revComp(sequences[i]))
    tempfile.close()
    muscleout= os.popen("./muscle " + quiet +"< temp.fas").read() #store muscle output
    # parse muscle output
    sequenceMatches = re.findall("\n[actgACTG\-\n]+",muscleout)
    nameMatches = re.findall(">.+",muscleout)
    tempReads = [] # store the two reads, one for forward direction the other for reverse
    tempNames = [] # store the two names for the two reads, one for forward and the other for reverse
    # reference sequence is not necessarily the first sequence in output
    for j in range(len(nameMatches)):
        if nameMatches[j][1:] == referenceName:
            references.append(sequenceMatches[j].replace("\n",""))
            currentReference = sequenceMatches[j]
        else:
            tempNames.append(nameMatches[j][1:])
            tempReads.append(sequenceMatches[j].replace("\n",""))
    # score the two reads and decide which is better
    score = 0 # for matches in the first read increment score, for matches in the second read decrement score
    for k in range(len(currentReference)):
        if currentReference[k].upper() == tempReads[0].upper():
            score = score + 1
        else:
            score = score - 1
        if currentReference[k].upper() == tempReads[1].upper():
            score = score - 1
        else:
            score = score + 1
    if score >= 0:
        sequenceNames[i] = tempNames[0]
        reads.append(tempReads[0])   
    else:
        sequenceNames[i] = tempNames[1]
        reads.append(tempReads[1])  
    i = i + 1
os.remove("temp.fas")

# create a composite reference sequence
# reference alignments cannot differ in the order of letters [agctAGCT], but they can differ by where gaps are inserted
# whenever a reference alignment has a gap that the others don't have, introduce that same gap into both the reference and the corresponding read
compositeReference = references[0] # the composite reference that we'll create by resolving differences between reference alignments
i = 1
while i < len(references):
    #compare references
    j = 0
    while j < max(len(compositeReference), len(references[i])):
        maxLength = max(len(compositeReference), len(references[i]))
        if j >= len(compositeReference):
            # other reference is longer, add gaps to end of compositeReference
            compositeReference = compositeReference + "-"
            reads[0] = reads[0] + "-"        
        elif j >= len(references[i]):
            # other reference is shorter, add gaps to end of references[i]            
            references[i] = references[i] + "-"
            reads[i] = reads[i] + "-"
        #print str(j)+","+str(len(compositeReference))+","+str(len(references[i]))
        compositeReference[j]
        references[i][j]
        if compositeReference[j] == "-" and not references[i][j] == "-":
            # introduce new gap into reads[i] at position j
            reads[i] = reads[i][0:j] + "-" + reads[i][j:]
            # introduce new gap into reference[i] at position j
            references[i] = references[i][:j] + "-" + references[i][j:]
            # add gap to end of compositeReference and corresponding read
            if len(references[i]) > len(compositeReference):
                compositeReference = compositeReference + "-"
                reads[0] = reads[0] + "-"
        elif not compositeReference[j] == "-" and references[i][j] == "-":
            # introduce new gap into compositeReference at position j
            compositeReference = compositeReference[0:j] + "-" + compositeReference[j:]
            # add gap to end of current read and reference
            if len(compositeReference) > len(references[i]):
                references[i] = references[i] + "-"
                reads[i] = reads[i] + "-"
            for k in range(i-2):
                # introduce new gap into all reads before read i
                reads[k+1] = reads[k+1][0:j] + "-" + reads[k+1][j:]
                # introduce new gap into reference[k] at position j
                references[k+1] = references[k+1][:j] + "-" + references[k+1][j:]
        j = j + 1
    i = i + 1
#print compositeReference
#print reads[0]
#print reads[1]
#print reads[2]
referenceSequence = compositeReference + "-"*(maxLength-len(compositeReference))
# add additional gaps in case any sequnece is too large
for i in range(len(reads)):
    reads[i] = reads[i] + "-"*(maxLength-len(reads[i]))
# compute coverage and flags 
flagString = ""
coverageString = ""
readSequences =[""]*len(reads) 
for i in range(len(compositeReference)):
    flagCount = 0
    coverageCount = 0;
    for j in range(len(reads)):
        #print "comparing " +reads[j][i].upper() + " and " + compositeReference[i].upper()
        if reads[j][i].upper() == compositeReference[i].upper() and not compositeReference[i] == "-":
            coverageCount = coverageCount + 1;
            readSequences[j] = readSequences[j] + reads[j][i].upper()
        else:
            if not compositeReference[i].upper() == "-":
                flagCount = flagCount + 1
                readSequences[j] = readSequences[j] + reads[j][i].lower()
            else:
                readSequences[j] = readSequences[j] + reads[j][i].lower()
    if coverageCount > 0:
        coverageString = coverageString + str(coverageCount)
    else:
        coverageString = coverageString + " " # append blank space instead of 0s to improve readability
    if flagCount > 0:
        flagString = flagString + str(1)
    else:
        flagString = flagString + " "
# print output
spacing = maxNameLength + 6
print ('{0: <' + str(spacing) + '}').format("Coverage") + "\t" + coverageString
print ('{0: <' + str(spacing) + '}').format("Flags") + "\t" + flagString
print ('{0: <' + str(spacing) + '}').format(referenceName) + "\t" + compositeReference
for i in range(len(sequenceNames)-1):
    print ('{0: <' + str(spacing) + '}').format(sequenceNames[i+1]) +"\t" + readSequences[i]

















