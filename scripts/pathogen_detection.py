import getopt
import multiprocessing
import os
import sys
import subprocess
import math
from colorama import Fore
from os import listdir
from os.path import isfile, join
from bioservices import *

__author__ = 'basir'
temp_dir = 'tmp'
pathogen_dir = 'pathogens'
ass_dir = 'assembly'
colorful = True


def cmd_exists(cmd):
    return subprocess.call('type ' + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0


def write_command(message):
    if colorful:
        print Fore.GREEN + '>>>>> ' + message + Fore.RESET
    else:
        print '>>>>> Command: ' + message


def write_warning(message):
    if colorful:
        print Fore.YELLOW + message + Fore.RESET
    else:
        print '>>>> Warning: ' + message


def write_error(message):
    if colorful:
        print Fore.RED + message + Fore.RESET
    else:
        print '>>>> Error: ' + message


def write_info(message, mode='m'):
    if mode == 'M':
        if colorful:
            print '>>> ' + message
        else:
            print '>>>> Info: ' + message
    else:
        if colorful:
            print Fore.YELLOW + ' ' * 3 + message + Fore.RESET
        else:
            print ' ' * 3 + '>>>> Info: ' + message


def process_parameters():
    '''
    Pathogen Detection Tool 0.0.1
    You should at least provide 5 input arguments for this tool:
    1- Provide the reference genomes files in FASTA format with value -r arguments
    2- Provide the pathogen database file in FASTA format with value -p argument.
    3- Provide the forward reads file in FASTQ format with value -f argument.
    4- Provide the reverse reads file in FASTQ format with value -b argument.
    5- Provide the output folder for the final report of tool with -o argument.
    ----------------------------------------------------------------------------
    Note that this tool needs bwa 0.6.2 and samtools 0.1.19 in order to run.
    '''
    refs = []
    pathogen_db = []
    fr = ''
    rr = ''
    output = ''
    clean= False
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hcr:p:f:b:o:',
                                   ['help', 'clean', 'reference=', 'pathogens=', 'forward_reads=', 'reverse_reads=', 'output='])
    except getopt.error, msg:
        write_info(msg)
        write_info('for help use --help')
        sys.exit(2)
        
        # process options
    for option, value in opts:
        if option in ('-h', '--help'):
            print process_parameters.__doc__
            sys.exit(0)
        if option in ('-c', '--clean'):
            clean=True
    
        if option in ('-r', '--reference'):
            if os.path.isfile(value) and os.stat(value)[6] != 0:
                refs.append(value)
            else:
                write_error('The reference genome file {0} you provided does not exist!'.format(value))
                sys.exit(1)
        if option in ('-f', '--forward_reads'):
            if os.path.isfile(value):
                fr = value
            else:
                write_error('The forward reads file {0} you provided does not exist!'.format(value))
                sys.exit(1)
        if option in ('-b', '--reverse_reads'):
            if os.path.isfile(value):
                rr = value
            else:
                write_error('The reverse reads file {0} you provided does not exist!'.format(value))
                sys.exit(1)
        if option in ('-p', '--pathogens'):
            if os.path.isfile(value):
                pathogen_db.append(value)
            else:
                write_error('The pathogens genome file {0} you provided does not exist!'.format(value))
                sys.exit(1)
        if option in ('-o', '--output'):
            if os.path.isdir(value):
                output = value
            else:
                output = value
                write_info('The output directory {0} you provided does not exist. Creating it...'.format(value))
                os.mkdir(output)
            if not os.path.isdir(output + '/' + temp_dir):
                os.mkdir(output + '/' + temp_dir)
            if not os.path.isdir(output + '/' + temp_dir + '/' + pathogen_dir):
                os.mkdir(output + '/' + temp_dir + '/' + pathogen_dir)
            if not os.path.isdir(output + '/' + temp_dir + '/' + ass_dir):
                os.mkdir(output + '/' + temp_dir + '/' + ass_dir)

    for ref in refs:
        reference_folder = output + '/' + temp_dir + '/' + os.path.splitext(os.path.basename(ref))[0]
        if not os.path.isdir(reference_folder):
            os.mkdir(reference_folder)

    return [refs, fr, rr, pathogen_db, output, clean]


def not_already_generated(file_name):
    return not os.path.isfile(file_name) or os.stat(file_name)[6] == 0


def index_reference_genomes(refs, path_db):
    # indexing the reference genomes
    write_info('Indexing the reference genomes...', 'M')
    for genome in refs:
        if not_already_generated(genome + '.bwt'):
            write_info('Indexing the reference genome {0} ...'.format(genome))
            command = 'bwa index -a bwtsw ' + genome
            execute_command(command)
        else:
            write_info('Reference genome file {0} is already indexed.'.format(genome))
    # Indexing the pathogens
    write_info('Indexing the pathogen database...', 'M')
    if not_already_generated(path_db[0] + '.bwt'):
        write_info('Indexing the pathogen database {0} ...'.format(path_db))
        command = 'bwa index -a is ' + path_db[0]
        execute_command(command)
    else:
        write_info('Pathogen database file {0} you provided is already indexed!'.format(path_db))


def check_prerequisites(executables):
    for executable in executables:
        if not cmd_exists(executable):
            write_error('{0} is not accessible! exiting...'.format(executable))
            write_error('This tool need bwa, idba_ud and samtools to be installed and be accessible in command line.\n'
                        'If you have already installed them, please make sure to add their executable path to the '
                        'PATH environment variable')
            sys.exit(1)
    return True


def filename_wo_ext(filename):
    return os.path.splitext(os.path.basename(filename))[0]


def execute_command(command):
    write_command(command)
    subprocess.call(command, shell=True)

def getAllMappings(fields):
    #XA Flag: XA:Z:ref1,+1,126M,0;ref2,+1,126M,0;
    #mapping tuple config (chrm,genomeStart,cigar)
    mappings=[]
    #get reported read
    mappings.append((fields[2],int(fields[3])-1,fields[5]))
    #get multiple hits FOR BWA ONLY
    for i in range(10,len(fields)):
        if "XA:Z" in fields[i]:
            parseString=fields[i].split(":")[-1]
            alns = parseString.split(";")
            for aln in alns:
                 tokens = aln.split(",")
                 if not len(tokens) <= 1:
                     mappings.append((tokens[0],int(tokens[3]),tokens[2]))
    return mappings
    
    
def addToCoverageMap(mapping,alignmentCounts,depthTables):
    #grab name, starting pos, and cigar
    chrm = mapping[0]
    genomeStart = int(mapping[1])
    cigarString = mapping[2]
    alignmentRepresentation = ""
    integerString=""
    #Generate string to aid with depth counts format: m = match, _ = deletion in reference
    for char in cigarString:
        if char == "M":
            alignmentRepresentation+="m"*int(integerString)
            integerString = ""
        elif char == "I" or char == "N" or char == "S" or char == "H" or char == "P": 
            integerString = ""
        elif char == "D": 
            alignmentRepresentation+="_"*int(integerString)
            integerString = ""
        else:
            integerString+=char
    #append to depth position on the correct chromosome        
    for i in range(len(alignmentRepresentation)):
        if alignmentRepresentation[i]=="m":
	    depthTables[chrm][genomeStart+i]+=1
    #increment Counts    
    alignmentCounts[chrm]+=1
	
    
def dumpFastqFiles(readsToFastq,fqReads):
    #dumps tuples into a fastq file
    readCounts = 0
    fqHandle = open(fqReads,"w")
    for key in readsToFastq:
        if len(readsToFastq[key])==2:
            fqHandle.write("@"+key+"_1\n"+readsToFastq[key][0][0]+"\n+\n"+readsToFastq[key][0][1]+"\n")
	    fqHandle.write("@"+key+"_2\n"+readsToFastq[key][1][0]+"\n+\n"+readsToFastq[key][1][1]+"\n")
	    readCounts+=2
	else:
	    fqHandle.write("@"+key+"\n"+readsToFastq[key][0][0]+"\n+\n"+readsToFastq[key][0][1]+"\n")
	    readCounts+=1
    fqHandle.close()
    return readCounts
	    
	    
	    
def calcStats(depthTables,alignmentCounts,isChrom,ResFile):
    # for chromosomes, Combine all contigs for all chromosomes
    resHandle = open(ResFile,"a")
    resHandle.write("Name\tLength\tAverage Depth\tSD of Depth\tCoverage\tTotal Mappings\n")
    totalReadsThatMapped=0
    if isChrom:
        tempTableCounts = {}
        tempTableDepth = {}
        for key in alignmentCounts:
	    chrom=key.split("_")[0]
	    if not chrom in tempTableCounts:
		tempTableCounts[chrom]=0
		tempTableDepth[chrom]=[]
            tempTableCounts[chrom] += alignmentCounts[key]
	    tempTableDepth[chrom] += depthTables[key]

	alignmentCounts =  tempTableCounts
	depthTables = tempTableDepth
	
    # Get stats
    for key in depthTables:
        if alignmentCounts[key] == 0:
            continue
        #avg depth
        avgDepth = sum(depthTables[key])/ float(len(depthTables[key]))
        #coverage
        coveredBp = 0
        for coverageVal in depthTables[key]:
	    if not coverageVal == 0:
	        coveredBp+=1
	coverage = coveredBp/float(len(depthTables[key]))
	#sd
	sd = 0.0
	for coverageVal in depthTables[key]:
	    sd+=(avgDepth-coverageVal)**2
	sd/=len(depthTables[key])    
	sd = math.sqrt(sd)
	#append counts total
	totalReadsThatMapped+=alignmentCounts[key]
	
	resHandle.write( key+"\t"+str(len(depthTables[key]))+"\t"+str(avgDepth)+"\t"+str(sd)+"\t"+str(coverage)+"\t"+str(alignmentCounts[key])+"\n")
    resHandle.write("All Mappings: "+str(totalReadsThatMapped)+"\n")
    resHandle.close()
  
def processSamFile(binaryMask,samFileName,fastqPrefix,isChrom,ResFile):

    
    #SEQuence HEADER example:@SQ	SN:ref1	LN:253 
     
    #SAM FIELDS
    #1 Qname
    #2 FLAG
    #3 RNAME
    #4 POS
    #5 MAPQ
    #6 CIGAR
    #7 RNEXT
    #8 PNEXT
    #9 TLEN
    #10 SEQ
    #11 QUAL
     
    #XA Flag: XA:Z:ref1,+1,126M,0;ref2,+1,126M,0;
    
    #keeps track of depth per bp
    depthTables={}
    #counts reads that hit each reference
    alignmentCounts={}
    #tuples of reads to add into fastq of unmapped reads
    readsToFastq = {}
    samFile = open(samFileName,"r")
    fqReads = fastqPrefix+"_merged.fastq"
    
    
    for line in samFile:
        fields = line.split("\t")
        
        #New Chrome/Pathogen create depth table and count table for it
        if ("@SQ" in line):
            depthTables[fields[1].split(":")[1]] = [0]*int(fields[2].split(":")[1])
            alignmentCounts[fields[1].split(":")[1]] = 0 
            continue
	#Header we do not care about   
        if line[0] == "@": 
            continue
        #if alignment line doesn't have all 11 fields discard the line    
        if len(fields) < 11:
            continue    
        #Is the mask we sent in hit?
        if int(fields[1]) & binaryMask > 0:
            if not fields[0] in readsToFastq:
	        #read did not meet our criterea....
                readsToFastq[fields[0]] = []
            readsToFastq[fields[0]].append((fields[9].rstrip(),fields[10].rstrip()))
        else:
	    #read mapped according to our criteria
            mappings = getAllMappings(fields)
            #Get all reported mappings!
            for mapping in mappings:
                addToCoverageMap(mapping,alignmentCounts,depthTables)
                
    #gen fastq file           
    readCounts= dumpFastqFiles(readsToFastq,fqReads)
    #generate stats of reads that did map
    calcStats(depthTables,alignmentCounts,isChrom,ResFile)
    fileHandle=open(ResFile,"a")
    fileHandle.write("Reads that did not map: "+str(readCounts))
    fileHandle.close()
    return fqReads    

def align_reads_to_references(genomes, fr, rr, output_dir, clean):
	#noc = multiprocessing.cpu_count()
	noc = 20
	for genome in genomes:
		genome_dir = output_dir + '/' + temp_dir + '/' + filename_wo_ext(genome)
		if not os.path.isdir(genome_dir):
			os.mkdir(genome_dir)
		write_info('Aligning the reads to the reference genome {0} ... '.format(genome))
		forward_read_filename = genome_dir + '/' + filename_wo_ext(fr)
		sam_file = forward_read_filename + '.sam'
		bam_file = forward_read_filename + '.bam'
		cov_file = forward_read_filename + '.cov'
		sai1_file = forward_read_filename + '1.sai'
		sai2_file = forward_read_filename + '2.sai'
	        execute_command("bwa aln -t " + str(noc) + " " + genome + " " + fr + " > " + sai1_file)
		execute_command("bwa aln -t " + str(noc) + " " + genome + " " + rr + " > " + sai2_file)
		execute_command("bwa sampe " + genome + " " + sai1_file + " " + sai2_file + " " + fr + " " + rr + " > " + sam_file)
		os.remove(sai1_file)
		os.remove(sai2_file)
                readPrefix = forward_read_filename + '_' + filename_wo_ext(genome) 
	ret =  processSamFile(12,sam_file,readPrefix,True,(output_dir+"/genomeResults.txt"))
	if clean:
	    os.remove(sam_file)
	return ret


def align_reads_to_pathogens(genomes, reads, output_dir, clean):
    #noc = multiprocessing.cpu_count()
    print genomes
    print reads
    print output_dir
    noc = 20
    for genome in genomes:
        genome_dir = output_dir + '/' + temp_dir + '/pathogen'
        if not os.path.isdir(genome_dir):
            os.mkdir(genome_dir)
        write_info('Aligning the reads to the reference genome {0} ... '.format(genome))
        readPrefix = genome_dir + '/pathogen' 
        sam_file = readPrefix + '.sam'
        sai_file = readPrefix + '_merged.sai'
        
       
       
        execute_command("bwa aln -t {0} {1} {2} > {3}".format(noc,genome,reads, sai_file))
        execute_command("bwa samse -n 1000 {0} {1} {2} > {3}".format(genome,sai_file,reads,sam_file))
        os.remove(sai_file)
    ret = processSamFile(4,sam_file,readPrefix,False,(output_dir+"/pathogenResults.txt"))
    if clean:
        os.remove(sam_file)
        os.remove(ret)
        os.remove(reads)



def index_pathogen_d(output_dir):
    viral_directory = output_dir + '/' + temp_dir + '/' + pathogen_dir
    for directory in os.walk(viral_directory).next()[1]:
        if directory.startswith('.'):
            continue
        files = [f for f in listdir(viral_directory + '/' + directory) if
                 isfile(join(viral_directory + '/' + directory, f))]
        for fasta_file in files:
            if fasta_file.endswith('.fa'):
                full = viral_directory + '/' + directory + '/' + fasta_file
                if os.path.exists(full + '.bwt'):
                    write_info(
                        'Pathogen {0} from family {1} has already been indexed!'.format(filename_wo_ext(fasta_file),
                                                                                        directory))
                    continue;
                write_info('Indexing the pathogen {0} from family {1}'.format(filename_wo_ext(fasta_file), directory))
                subprocess.call(['bwa', 'index', full, '-a', 'is'])


def assemble_remaining_reads(filteredReads, output_dir):
    write_info('Assembling the remaining reads to pathogens database...', 'M')
    merged_fasta = output_dir + '/' + temp_dir + '/' + filename_wo_ext(filteredReads) + '.fa'
    command = 'bin/fq2fa --filter ' + filteredReads + ' ' + merged_fasta
    write_command(command)
    os.system(command)
    assembly_output = output_dir + '/' + ass_dir + '/' + filename_wo_ext(filteredReads)
    
    if not os.path.isdir(output_dir + '/' + ass_dir):
        os.mkdir(output_dir + '/' + ass_dir)
    if not os.path.isdir(assembly_output):
        os.mkdir(assembly_output)
    execute_command('bin/idba_ud --num_threads 12 -o ' + assembly_output + ' -r ' + merged_fasta)
    if os.path.isfile(assembly_output+"/scaffold.fa"):
        os.system("mv "+assembly_output+"/scaffold.fa "+output_dir)
        os.system("gzip "+output_dir+"/scaffold.fa")
    elif os.path.isfile(assembly_output+"/contig.fa"):
        write_warning("Couldn't find assembled scaffolds reporting contigs instead")
        os.system("mv "+assembly_output+"/contig.fa "+output_dir)
        os.system("gzip "+output_dir+"/contig.fa")
    else:
        write_warning("Warning: Couldnt find contigs or scaffolds from assembly")
    os.system("rm -rf "+output_dir + '/' + ass_dir)    
     


def main():
    check_prerequisites(['bwa', 'samtools', 'bin/idbaud/bin/idba_ud'])
    [refs, fr, rr, path_db, output_dir, clean] = process_parameters()
    index_reference_genomes(refs, path_db)
    filteredReads = align_reads_to_references(refs, fr, rr, output_dir, clean)
    assemble_remaining_reads(filteredReads, output_dir)
    align_reads_to_pathogens(path_db, filteredReads, output_dir, clean)



if __name__ == '__main__':
    main()


