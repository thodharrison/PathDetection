from bioservices import *
from Bio.Blast import NCBIWWW

#sequence = UniProt().searchUniProtId("P43403", "fasta")
#print type(sequence)

#print sequence
sequence = ""
with open ("results/assembly/1_ATCACG_L006_R1_001_felCat5_1_pathogens/contig.fa", "r") as myfile:
    sequence=myfile.read()

result_handle = NCBIWWW.qblast("blastn", "nt", sequence,format_type="Text")
print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"	
print len(sequence)
save_file = open("my_blast.xml", "w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()

#s = NCBIblast(verbose=True)

#jobid = s.run(program="blastn", sequence=sequence, 
#              stype="dna", database="", 
#              email="b.shariat@gmail.com")



#print(s.getResult(jobid, "out")[0:1200])


