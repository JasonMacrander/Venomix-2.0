#!/usr/bin/python
"""
Programing Suite: Venomix(xxx) - Main

This file retrieves information from the BLAST output and the 
ToxProt information file ().
The script first takes the UniProt table (uniref_50), but you can upload
your own tabular dataset relative to the query sequences.
together by linking "funcitonal groups" with Transcripts.

Version 0.7 (Jason Macrander)
"""
import re
import sys
sys.path.insert(0,'support_files')
import numpy
from Function import RSEMtoDict
from Function import Combine_Files
from Function import Convert_Blast
from Function import Tox_Groups
from Function import Genbank_Summary
from Function import Toxinfiles
from Function import genbankinfo
from Function import expressioninfo
from subprocess import call
from Bio.Align.Applications import MafftCommandline
from io import StringIO
from Bio import AlignIO 
import getopt, os, subprocess, sys
from Alignment import Align
import shutil
import getopt
from distutils.spawn import find_executable

try:
    from Bio import SeqIO
except:
    print("This program requires the Biopython library")
    sys.exit(0)

from Venomixxxxx import BEGIN

#Input files provided by the users

#<<<<<<< HEAD
Tox_Cluster = "../data/ToxProt/uniref_50.txt" #Standard
gbkFile = "../data/ToxProt/sequence.gp" #Standard

#=======
Blast_out_name=''
RSEM_fileinput=''
Trans_input=''
mulexp=''
#1myopts, args = getopt.getopt(sys.argv[1:],"b:r:f:x")
myopts, args = getopt.getopt(sys.argv[1:],"b:d:f:r:x")
for o, a in myopts:
 if o == '-b':
  Blast_out_name=("{}".format(a))
 elif o == '-d':
  RSEM_fileinput=("{}".format(a))
 elif o == '-f':
  Trans_input=("{}".format(a))
 elif o == '-r':
  mulexp = a.split(',')
 elif o == '-x':
  BEGIN()

Full_Trans=open(Trans_input, 'r')
a=os.getcwd()
#>>>>>>> 2e363297f9570ce8072da323499b38a83711ef93

#Intermediate files for your own interest/use
TPMFile = open("../../TPM.fasta", 'w')
ToxinLike = ("../../Toxin_like.txt")
test = 0
#BEGIN(test)

"""
This function retrieves information from the RSEM file and Fasta Assembly
and writes to a file based on certain parameters.
"""

if RSEM_fileinput != '':
	RSEM_file=open(RSEM_fileinput, 'r')
	RSEM_Dict = {}
	RSEM_Dict = RSEMtoDict(RSEM_file)
	print("RSEM Dict Loaded!")
	RSEM_Seq_Dict = Combine_Files(RSEM_Dict,Full_Trans)
	print("Seq Dict Loaded!")
	count = 0
	for key in RSEM_Seq_Dict:
		TPM = float(RSEM_Seq_Dict.get(key).get('TPM'))
		if (float(RSEM_Seq_Dict.get(key).get('TPM')) > 1.0):
			TPMFile.write(">%s\n%s\n" % (key , RSEM_Seq_Dict.get(key).get('Sequence')))
			TPM_Trans = (">%s\n%s\n" % (key , RSEM_Seq_Dict.get(key).get('Sequence')))	
			count += 1
		else:
			continue
	#print "no errors?"

	print("Finished RSEM Step")

print("Inputs look good!")


"""
Begin the step which uses the Convert_Blast function to creat two dictionaries,
one with all of the ToxProt IDs as keys (TPM_Trans IDs are the eleements)and the 
other with all the TPM_Trans IDs as keys (ToxProt IDs are the eleements)
"""
print("Beginning Convert_Blast Step")

QueryD = Convert_Blast(Blast_out_name)

b = Genbank_Summary(gbkFile, Tox_Cluster)

with open("../../genbank.info", 'w') as f:
	for item in b:
		f.write(str(item)+'\n')
f.close()

#Input the Transcriptomes file
fasta_sequences = list(SeqIO.parse(Trans_input, 'fasta'))

Groups = {}
(Groups, Genbank_Summary, Genbank_Full) = Tox_Groups(Tox_Cluster, QueryD, Blast_out_name, gbkFile)

# Make directories for temp files involving subsequent work being done
seq_list = "../data/Inputs/temp/list"
if not os.path.exists(seq_list):
	os.makedirs(seq_list, mode=0o777)

	
seq_out = "../data/Inputs/temp/out"
if not os.path.exists(seq_out):
	os.makedirs(seq_out, mode=0o777)

seq_error="../data/Inputs/temp/Error"
if not os.path.exists(seq_error):
	os.makedirs(seq_error, mode=0o777)
	
#Creating folders naming every toxin
#<<<<<<< HEAD

Toxinfiles(Groups, ToxinLike, fasta_sequences)

"""
#This step will use the program Transdecoder to on target sequences
"""
#>>>>>>> 2e363297f9570ce8072da323499b38a83711ef93
Transdir="../code/Transdecoder"
if not os.path.exists(Transdir):
	os.makedirs(Transdir, mode=0o777)


print ("Transdecoder step started:")
os.chdir(Transdir)
for filename in os.listdir("../../data/Inputs/temp/out/"):
	if filename.endswith(".txt"):
		#print filename
		name = os.path.join("../../data/Inputs/temp/out/", filename)
		with open("../../data/Inputs/temp/Error/{}".format(filename), 'w') as f:
			subprocess.call('../../Programs/Trans/TransDecoder.LongOrfs -t %s' %name, shell=True, stdout=open(os.devnull,'wb'), stderr=f)
			subprocess.call('../../Programs/Trans/TransDecoder.Predict -t %s' %name, shell=True, stdout=open(os.devnull,'wb'), stderr=f)
print("Transdecoder step completed")
os.chdir(a)
#Input uniprot database
fastauniprot_sequences = list(SeqIO.parse("../data/ToxProt/uniprot_ID.fasta",'fasta'))

#Create a finaloutput directory for further use.
outputdir="../FinalOutput"
if not os.path.exists(outputdir):
	os.makedirs(outputdir, mode=0o777)

#Run the Alignment script on all the output directories of Transdecoder
#Run genbankinfo and expressioninfo on the output files.
Align(outputdir)
for file in os.listdir(outputdir):
	name = os.path.join(outputdir, file)
	for filename in os.listdir(name):
		if filename == "uniprot.fasta":
			name1 = os.path.join(name,filename)
			genbankinfo(name1, name)
		if filename == "sequence.fasta":
			name1 = os.path.join(name, filename)
			for rsemfiles in mulexp:
				if rsemfiles != '':
					expressioninfo(name1, name, rsemfiles)

#Delete the transdecoder directories created earlier
shutil.rmtree(Transdir)
shutil.rmtree("../data/Inputs")
for dir in os.listdir(outputdir):
	path1 = os.path.join(outputdir, dir)
	if os.path.isdir(path1) and len(os.listdir(path1)) == 0:
		os.rmdir(path1)

subprocess.call('rm support_files/*.pyc', shell=True)

print("Venomix has completed")