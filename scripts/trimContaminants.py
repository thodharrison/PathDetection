import sys , os
# usage <input csv> <Output folder>
if not len(sys.argv) == 3:
    print "usage: \npython scripts/trimContaminants.py <input.csv> <Output Folder>"
    sys.exit()
argFile = open(sys.argv[1].rstrip(),"r")
outdir = sys.argv[2].rstrip()
fileTups=[]
for line in argFile:
    fileVec = line.split(",")
    if len(fileVec) == 2:
        print "Detected Single Ended reads, The pipeline is only written for paired end reads:" 
        sys.exit()
    elif len(fileVec) == 3:
        if not fileVec[0] == "":
            fileTups.append((fileVec[0].rstrip(),fileVec[1].rstrip(),fileVec[2].rstrip()))
            
    if not os.path.exists(outdir+"/pairedReads"):
        os.makedirs(outdir+"/pairedReads")
    if not os.path.exists(outdir+"/unpairedReads"):
        os.makedirs(outdir+"/unpairedReads")         
    
    for tup in fileTups:
        nameFr = tup[0].split(".")[0].split("/")[-1]
        nameRr = tup[1].split(".")[0].split("/")[-1]
        os.system("java -jar bin/Trimmomatic-0.33/trimmomatic-0.33.jar PE "+tup[0]+" "+tup[1]+" "+outdir+"/pairedReads/"+nameFr+"_trimmed_fq.gz "+outdir+"/unpairedReads/"+nameFr+"_trimmed_fq.gz "+outdir+"/pairedReads/"+nameRr+"_trimmed_fq.gz "+outdir+"/unpairedReads/"+nameRr+"_trimmed_fq.gz ILLUMINACLIP:"+tup[2]+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >& "+outdir+"/"+nameFr+".txt ")
        #print "java -jar bin/Trimmomatic-0.33/trimmomatic-0.33.jar PE "+tup[0]+" "+tup[1]+" "+outdir+"/pairedReads/"+nameFr+"_trimmed_fq.gz "+outdir+"/unpairedReads/"+nameFr+"_trimmed_fq.gz "+outdir+"/pairedReads/"+nameRr+"_trimmed_fq.gz "+outdir+"/unpairedReads/"+nameRr+"_trimmed_fq.gz ILLUMINACLIP:"+tup[2]+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >& "+outdir+"/nameFr.txt &"
        #print "######"

'''
readDir = sys.argv[1]
if not os.path.exists(directory):
    os.makedirs(directory)
names = os.listdir(sys.argv[1])
fqs=[]
for name in names:
    if ".fq"in name or ".fastq" in name:
        fqs.append(name)
names.sort()
print names        
'''         




