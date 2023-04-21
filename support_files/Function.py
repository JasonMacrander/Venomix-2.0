#!/usr/bin/python
"""
Programing Suite: Venomix(xxx) - Functions
Version 0.7 (Jason Macrander)
"""
import sys
import re
import numpy
from Bio import SeqIO
from Bio import GenBank
from Bio.SeqIO.FastaIO import SimpleFastaParser
from distutils.spawn import find_executable

"""
loop through each line of the RSEMFile, determine if it is a 
.genes or .isoforms file, put the items from that file into a dictionary.
"""
def RSEMtoDict(RSEMFile):
	Info_dic = {}
	LineCount = 0
	for Line in RSEMFile:
		LineCount += 1
		if LineCount == 1:
			continue
		if LineCount >= 2:	
			Line = Line.strip('\n')
			#If the file looks like a .genes output.			
			ElementList = (Line.split('\t'))
			RSEM_dic = {"ID": ElementList[1], 
			"Length": ElementList[2], 
			"E_length" : ElementList[3], 
			"E_count" : ElementList[4], 
			"TPM" : ElementList[5], 
			"FKPM" : ElementList[6]}
			if Line.count('\t') == 6:
				RSEM_dic["IsoPct"] = None
				RSEM_dic["Type"] = "gene"
			elif Line.count('\t') == 7:
				RSEM_dic["IsoPct"] = ElementList[7]
				RSEM_dic["Type"] = "isoform"
			else:
				print("Unknown file type.")
				return
			Info_dic[ElementList[0]] = RSEM_dic
	return Info_dic

"""
Comment about function
"""
def Combine_Files (RSEM,SeqFIle):
	with SeqFIle as handle:
		for values in SimpleFastaParser(handle):
			Fasta_ID = values[0].partition(" ")[0]
			if Fasta_ID == 'comp1000001_c0_seq1': 
				continue
			else:
				RSEM[Fasta_ID]["Sequence"] = values[1]
	return RSEM	

"""
Comment about function
"""
def Convert_Blast(Blast_out_name):
	Blast_Out = open(Blast_out_name, 'r')
	QueryD = {}
	for line in Blast_Out:
		split_line = line.split('\t') #split tabular blast output
		query, hit = split_line[:2] #retrieve query information from line
		focus = query.split("|")[1] #Isolate UniProt ID for downstream matching
		QueryD[focus] = QueryD.get(focus,[]) + [hit]
	return (QueryD)

"""
This file retrieves information from the BLAST output and the 
ToxProt information file (). The script first takes the UniProt 
table (uniref_50), but you can upload your own tabular dataset 
relative to the query sequences. together by linking 
"funcitonal groups" with Transcripts.
"""
def Tox_Groups (Tox_Cluster, QueryD, Blast_out_name, gbkFile):
	gbk = open(gbkFile, 'r')
	Blast_Out = open(Blast_out_name, 'r')
	Tox_File = open(Tox_Cluster, 'r')
	QueryKeyList = QueryD.keys()
	Groups = {}
	Genbank_Summary = {}
	Genbank_Full = {}
	for record in SeqIO.parse(gbk, "genbank"):
		gb_ID = record.id.split(".")[0]
		Genbank_Full[gb_ID] = record
		temp_Desc = (record.description)
		temp_Desc = gb_ID + ": " + temp_Desc[14:]
		temp_Desc = temp_Desc.split(";")[0]
		Description = [' '.join(temp_Desc.split())]		
		Genbank_Summary[gb_ID] = Genbank_Summary.get(gb_ID, []) + list(set(Description)) 
											
	LineCount = 0
	Gene_group = []
	for line in Tox_File:
		LineCount = LineCount + 1
		keyCount = 0
		keyCount1 = 0
		if LineCount == 1:
			print("Skipping Header")
		if LineCount >= 2:
			ACCESSION = line.split("\t")[4].split("; ")
			for ID in ACCESSION:
				ACC = ID.split("; ")
			GeneC = line.split("\t")[2].split(": ")
			Gene_group = ''.join(GeneC[1:])
			Gene_group = Gene_group.replace(" (Fragments)", "")
			Gene_group = Gene_group.replace(" (Fragment)", "")
			numbers = line.split("\t")[3]
			Gene_key = "%s(%s)" % (Gene_group,numbers)
			if int(numbers) > 1:
				keyCount == keyCount + 1
			for Akey in ACC:
				if Akey in QueryKeyList:
					keyCount1 = keyCount1 + 1
					Groups[Gene_key] = Groups.get(Gene_key,[]) + list(set(QueryD[Akey]))
	print("Printing to Summary File")
	print("Printing to Full File")
	print("Toxin gene groups: ", len(Groups.keys())) 
	Blast_Out.close()
	Tox_File.close()
	return (Groups, Genbank_Summary, Genbank_Full)

"""
This function retrieves information from the ToxProt.txt file 
and writes to a dictionary, grouping the file based on the 
left line indicators.
"""
def Genbank_Summary(gbkFile, Tox_Cluster):
	gbk = open(gbkFile, 'r')
	LineCount = 0
	Tox=open(Tox_Cluster, 'r')
	a = []
	for line in Tox:
		LineCount = LineCount + 1
		if LineCount == 1:
			print("Skipping Header")
		if LineCount >= 2:	
			for record in SeqIO.parse(gbk, "genbank"):
				a.append(record)
	#print("COMPLETED")
	return a

"""
This function retrieves information from the Tox_Groups output and the 
to categorize candidate toxin genes into a directory (date specific).
Fasta files for each cluster of toxin candidates are labeled according 
to their associative toxin groups.
"""

def Toxinfiles(groups, toxinlike, fasta):
	fasta_sequences=fasta
	ToxinLike = open(toxinlike, 'w')
	Groups=groups
	for item in Groups:
		for key in Groups.keys():
			ToxinLike.write("\n%s:\n%s" % (item, Groups[item]))	
		
		item1=item.replace("/", "@")
		item2=item1.replace("(", "_")
		item2=item2.replace(")", "")
		item2=item2.replace(" ", "_")
		wanted = set()	
		with open("../data/Inputs/temp/list/{}.txt".format(item2), 'w') as file:
			for line in Groups[item]:
				line.strip('[]')
				file.write("\n{}".format(line))			
				if line != "":
					wanted.add(line)
			
		frags=[]
		with open("../data/Inputs/temp/out/{}.txt".format(item2), 'w') as f:
			for seq in fasta_sequences:
				name = seq.id
				if name in wanted and len(seq.seq) > 0:
					wanted.remove(name)
					frags.append(seq)
			SeqIO.write(frags, f, "fasta")

'''
The genbank file is created using sequence.gp. 
'''
def genbankinfo(infile, outfile):
	f3=[]
	with open(infile, 'r') as f:
		for line in f:
			if line.startswith(">"):
				f3.append(line.strip('>').strip('\n'))
	with open("{}/genbankfile.txt".format(outfile), 'w') as f4:
		for item in f3:	 
			with open("../data/ToxProt/sequence.gp", 'r') as f1:
				for line in f1:
					if item in line:
						break
				for line in f1:
					f4.write(line)				
					if "//\n" in line:
						break		
						
						
def expressioninfo(infile, outfile, rsemfile):
	file1=open(rsemfile, 'r')
	infile=open(infile, 'r')
	a=[]
	i=0
	for line in infile:
		if line.startswith('>'):
			a.append(line.split('::')[1])
			#line1 = line.rstrip()
			#a.append((line.replace('>','')).rstrip())
	#print(a)
			#a.append(line.split('>')[0])
			#a.append(line.split('\n')[0])

	with open("{}/rsemfile_expressioninfo.txt".format(outfile), 'a') as f1: 
		for line in file1:
			if i == 0:
				f1.write(line)
			i += 1

		for item in a:
			with open(rsemfile, 'r') as file1:
				for line in file1:
					if item in line:
						RSEM_OUTPUT = rsemfile + "\t" + line
						f1.write(RSEM_OUTPUT)
