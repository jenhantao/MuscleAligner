#! /usr/bin/python

# Process: use muscle to do pairwise alignment of the reference to all of the short reads/sequences
# output: combine all of the aliignments and show coverage of the reference.  Then we can compute things like, where is there a gap in the coverage, and where do two different reads disagree?

# input file: fasta file, with first sequene being the reference
# second argument is whether or not to silence muscle, silenced by default, n to output everything
import os
import sys
import re
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
alignments = [];
while i < len(sequenceNames):
    #print ">" + referenceName + "\n" + referenceSequence + "\n>" + sequenceNames[i] + "\n" + sequences[i]
    tempfile = open('temp.fas','w+')
    tempfile.write(">" + referenceName + "\n" + referenceSequence + "\n>" + sequenceNames[i] + "\n" + sequences[i])
    tempfile.close()
    i = i + 1
    muscleout= os.popen("./muscle " + quiet +"< temp.fas").read() #store muscle output
    # parse muscle output
    sequenceMatches = re.findall("\n[actgACTG\-\n]+",muscleout)
    nameMatches = re.findall(">.+",muscleout)
    # reference sequence is not necessarily the first sequence in output
    if nameMatches[0][1:] == referenceName:
        alignments.append(sequenceMatches[1].replace("\n",""))
    else:
        alignments.append(sequenceMatches[0].replace("\n",""))        
os.remove("temp.fas")

# compute coverage and flags 
flagString = ""
coverageString = ""
alignPositions = [] # stores the position of a read aligned to the reference relative to the original reference sequence
for i in range(len(referenceSequence)):
    flagCount = 0
    coverageCount = 0; 
    for j in range(len(alignments)):
        #print "comparing " +alignments[j][i].upper() + " and " + referenceSequence[i].upper()
        if alignments[j][i].upper() == referenceSequence[i].upper():
            coverageCount = coverageCount + 1;
        else:
            if not referenceSequence[i].upper() == "-":
                flagCount = flagCount + 1
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
print ('{0: <' + str(spacing) + '}').format(referenceName) + "\t" + referenceSequence
for i in range(len(sequenceNames)-1):
    print ('{0: <' + str(spacing) + '}').format(sequenceNames[i+1]) +"\t" + alignments[i]





















