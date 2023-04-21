import re
import numpy
from subprocess import call
import sys
from Bio.Align.Applications import MafftCommandline
from io import StringIO 
from Bio import AlignIO 
import getopt, os, subprocess, sys

#Import the Biopython package
try:
    from Bio import SeqIO
except:
    print ("This program requires the Biopython library")
    sys.exit(0)

# Uniprot database is given as input.
fastauniprot_sequences = list(SeqIO.parse("../data/ToxProt/uniprot_ID.fasta",'fasta'))

'''This function performs the below steps
Takes Transdecoder files and runs blastp search, and creates fasta files. Filters the sequences from Uniref with the same ID and performs mafft alignment.
Creates a tree using Rpackage. Using signalp , we filter out the fasta files based on expression info.
'''
def Align(output):
    out=output
    for dir in os.listdir("../code/Transdecoder"):
        if dir.endswith("_dir") and not dir.startswith("."):    
            dir1=os.path.join("../code/Transdecoder", dir)
            dir2=dir.split(".")[0]
            if not os.path.exists("{}/{}".format(out, dir2)):
                os.mkdir("{}/{}".format(out, dir2), 0o777)
            dir3=os.path.join("{}/{}".format(out, dir2))
            for file in os.listdir(dir1):
    
                if file.endswith(".pep"):
                    file1=os.path.join(dir1, file)
                    print("Running reciprocal BLAST step")
                    subprocess.call('../Programs/blast/bin/blastp -query %s -db "../data/ToxProt/BLAST_uniref_50/BLAST_uniref_50" -outfmt 6 -out %s.txt -max_target_seqs 1 -evalue .0001' % (file1, file1), shell=True)
                    if os.path.getsize("{}.txt".format(file1)) > 0:
                        sequence=[]
                        newsequence=[]
                        seq1=[]
                        IDs=[]
                        frags1=[]           
                        Species=[]                  
                        with open("{}.txt".format(file1), 'r') as i:
                            for line in i:
                                sample=line.split('\t')[0]
                                sequence.append(sample)
                                newsequence.append(line.split('\t')[1])
                                for seq in newsequence:
                                    seq2=seq.split('|')[1]
                                    seq1.append(seq2)
                        
                        with open("../data/ToxProt/uniref_50.txt", "r") as uniref:
                            for line in uniref:
                                line1=line.split('\t')[0]
                                line2=line1.strip('UniRef50_')
                                if line2 in seq1:
                                    line3=line.split('\t')[4]
                                    line4=line3.split(';')
                                    line5=line.split('\t')[5]
                                    line6=line5.split(';')
                                    for line in line4:
                                        line7=line.strip(' ')           
                                        IDs.append(line7)
                                    
                                    for line in line6:
                                        line8=line.strip(' ')
                                        Species.append(line8)
                        with open("{}/Species.txt".format(dir3), 'w') as n: 
                            for item in Species:
                                n.write(item)
    
                        sequences=list(SeqIO.parse(file1,'fasta'))              
                        
                        with open("{}/uniprot.fasta".format(dir3), 'w') as f:
                            for seq in fastauniprot_sequences:
                                name = seq.id
                                if name in IDs and len(seq.seq) > 0:
                                    frags1.append(seq)
                            SeqIO.write(frags1, f, "fasta")         
                        f.close()
                        
                        for file in os.listdir("../data/Inputs/temp/out"):
                              if dir2 in file:
                                    with open("{}/nucleotide.fasta".format(dir3), 'w') as f:
                                          with open("../data/Inputs/temp/out/{}".format(file)) as fp:
                                               for line in fp:
                                                    f.write(line)
							
                        with open("{}/sequence.fasta".format(dir3), 'w') as f:
                            frags=[]
                            for seq in sequences:
                                name = seq.id
                                if name in sequence and len(seq.seq) > 0:
                                    sequence.remove(name)
                                    frags.append(seq)   
                            SeqIO.write(frags, f, "fasta")
                            
                        f.close()                                   
                        
                        with open("{}/uniprot.fasta".format(dir3), 'r') as f:
                            f1contents=f.read()
                        
                        with open("{}/sequence.fasta".format(dir3), 'r') as j:
                            f2contents=j.read()             
                        
                        with open("{}/final.fasta".format(dir3), 'w') as k:
                            k.write(f1contents+f2contents)
                            
                        with open("{}/final.fasta".format(dir3), 'r') as f:
                            counter=0
                            for line in f:
                                if '>' in line:
                                    counter=counter+1                       
                            if counter>1:                       
                                in_file = "{}/final.fasta".format(dir3)
                                mafft_cline = MafftCommandline("../Programs/mafft-linux64/mafft.bat", input=in_file)
                                print(mafft_cline)
                                stdout, stderr = mafft_cline() 
                                align = AlignIO.read(StringIO(stdout), "fasta")                 
                                with open("{}/aligned.fasta".format(dir3), "w") as handle:   
                                    handle.write(stdout)
                                handle.close()
											   
                        with open("../data/Inputs/temp/Error/Rerror.txt", 'w') as k:
                            subprocess.call("Rscript support_files/ape.r %s/aligned.fasta %s/finaltree.pdf" % (dir3, dir3), shell=True, stderr=k)
                            subprocess.call('../Programs/signalp -m %s/mature.fasta %s/final.fasta > %s/signalp.txt' %(dir3,dir3,dir3) , shell=True)
	