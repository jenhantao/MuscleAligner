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
reads = []; # array that stores the alignments of the reads
references = []; # array that stores the alignments of the reference sequence
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
        reads.append(sequenceMatches[1].replace("\n",""))
        references.append(sequenceMatches[0].replace("\n",""))
    else:
        reads.append(sequenceMatches[0].replace("\n","")) 
        references.append(sequenceMatches[1].replace("\n",""))  
    #print ">" + referenceName + "\n" + references[-1] + "\n>" + "read" + "\n" + reads[-1]         
os.remove("temp.fas")

# create a composite reference sequence
# reference alignments cannot differ in the order of letters [agctAGCT], but they can differ by where gaps are inserted
# whenever a reference alignment has a gap that the others don't have, introduce that same gap into both the reference and the corresponding read
compositeReference = references[0] # the composite reference that we'll create by resolving differences between reference alignments
i = 1
while i < len(references):
    #compare references
    j = 0
    while j < len(compositeReference):
        if compositeReference[j] == "-" and not references[i][j] == "-":
            # introduce new gap into reads[i] at position j
            reads[i] = reads[i][0:j] + "-" + reads[i][j:]
            # introduce new gap into reference[i] at position j
            references[i] = references[i][:j] + "-" + references[i][j:]
            # add gap to end of compositeReference and corresponding read
            compositeReference = compositeReference + "-"
            reads[0] = reads[0] + "-"
        elif not compositeReference[j] == "-" and references[i][j] == "-":
            # introduce new gap into compositeReference at position j
            compositeReference = compositeReference[0:j] + "-" + compositeReference[j:]
            # add gap to end of current read and reference
            references[i] = references[i] + "-"
            reads[i] = reads[i] + "-"
            for k in range(i-1):
                # introduce new gap into all reads before read i
                reads[k] = reads[k][0:j] + "-" + reads[k][j:]
                # introduce new gap into reference[k] at position j
                references[k] = references[k][:j] + "-" + references[k][j:]
        j = j + 1
    i = i + 1
referenceSequence = compositeReference
# compute coverage and flags 
flagString = ""
coverageString = ""
readSequences =[""]*len(reads) 
for i in range(len(compositeReference)):
    flagCount = 0
    coverageCount = 0;
    for j in range(len(reads)):
        #print "comparing " +reads[j][i].upper() + " and " + compositeReference[i].upper()
        if reads[j][i].upper() == compositeReference[i].upper():
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

















