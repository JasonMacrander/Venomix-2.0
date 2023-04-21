#!/usr/bin/python
"""
This file retrieves information from the ToxProt.txt file and writes to a dictionary, grouping the file based on the left line indicators.
"""
import Bio

gbkFile = open("../../data/Inputs/sequence.gp", 'r')

Fn=()
from Bio import SeqIO
from Bio import GenBank

Features = []

for record in SeqIO.parse(gbkFile, "genbank"):
	Fn = (len(record.features))
	Features.append(Fn)
print(max(Features))
print(min(Features))
print(sum(Features)/float(len(Features)))

print("Record %s has %i features" % (record.name,len(record.features)))
print(record.id)
print(record.name)
print(record.description)
#print record.comments

print(repr(record.seq))
gbkFile.close()
