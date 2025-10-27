#! /usr/bin/env python
# source: https://github.com/jessicarick/refbias_scripts/blob/master/sim_scripts/vcf2phylip.py3
# edited by Danny Jackson to accept vcf.gz files and phased genotypes

# python version 2.7.2+
# by Joana, script to convert vcf to phylip
# Some functions are recycled from the RAD python script written by Sam Wittwer

import sys, getpass, re, argparse

import gzip, io

def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, mode)

def parse_gt(sample_field):
    """Return tuple of allele indices as strings, or None if missing."""
    # sample_field like '0|1:35:99' or '0/1' or '.' or './.'
    gt = sample_field.split(":", 1)[0]
    if gt in (".", "./.", ".|."):
        return None
    gt = gt.replace("|", "/")
    parts = gt.split("/")
    if len(parts) == 1:                 # haploid -> duplicate
        parts = [parts[0], parts[0]]
    if not all(p.isdigit() for p in parts):
        return None
    return parts[0], parts[1]

# Define the neccessary functions:

# Function to extract header info: extractHeaderInfo
# takes a vcf file from line 1, goes through header, returns:
# 0[string]             header
# 1[list[string]]       individual IDs
# 2[int]                number of individuals in vcf file
# 3[int]                length of header
def extractHeaderInfo(input):
        linecounter = 0
        header = ""
        for line in input:
                linecounter += 1
                if re.match("^\#\#", line):
                        header+=line
                else:
                        header+=line
                        n_individuals = len(str(line.strip('\n')).split('\t'))-9
                        id_individuals = str(line.strip('\n')).split('\t')[9:]
                        break
        return header, id_individuals, n_individuals, linecounter

# Function to make all sample names of equal length (could be simplified)
# FillUp: accepts list, fills up list with optional character (default is " ") until they are of equal length
def fillUp(list, fill = " "):           # Checks length of entries in a list, fills up with specified fill string.
        returnlist = []
        for entry in list:
                if len(entry) < len(max(list, key=len)):
                        returnlist.append(entry+(len(max(list, key=len))-len(entry))*fill)
                else:
                        returnlist.append(entry)
        return returnlist

	
# Function to write the the lines in phylip format	
def writePhylipSequences(samplenames, sequences, outputdestination, writeref):
        if writeref:
                beginning = 0
                nsamples = str(len(samplenames))
        else:
                beginning = 1
                nsamples = str(len(samplenames)-1)
        nbases = str(len(sequences[0]))
        outputdestination.write(nsamples +"\t"+nbases+"\n")
        outstring = ""
        for i in range(beginning, len(samplenames)):
                outstring += samplenames[i]+"".join(sequences[i])
                outstring +="\n"
        outputdestination.write(outstring.strip("\n"))
	

# Parse the arguments provided
parser = argparse.ArgumentParser(description='Convert vcf file to phylip file format')

parser.add_argument('-i', '--input', dest='i', help="input file in vcf format [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="output file [required]", required=True)
parser.add_argument('-r', '--ref', action='store_true', help="if -r is specified, the reference sequence will be included in the phylip file)", default=False)
parser.add_argument('-f', '--fill', action='store_true', help="if -f is specified, all sites not in the vcf file will be printed as missing (N)", default=False)
parser.add_argument('-e', '--exclIndels', action='store_true', help="if -e is specified, indels are not printed (else replaced by N)", default=False)

args = parser.parse_args()
	
# Set the default values:
input = open_maybe_gzip(args.i, "rt")
output = open(args.o, "w")

writeref = args.ref
fill = args.fill	
noIndels = args.exclIndels
prev=100000000000000000000000000000000000

# How to convert vcf-style genotypes to single letters
AmbiguityMatrix = [["A","M","R","W","N"],["M","C","S","Y","N"],["R","S","G","K","N"],["W","Y","K","T","N"],["N","N","N","N","N","N"]]
CoordinatesDictionary = { "A":0 , "C":1, "G":2, "T":3, ".":4 }

def get_iupac(gt_tuple, altlist):
    a, b = gt_tuple  # strings "0"/"1"/...
    ai = int(a); bi = int(b)
    return AmbiguityMatrix[CoordinatesDictionary[altlist[ai]]][CoordinatesDictionary[altlist[bi]]]


# If -f and -e are specified 
if noIndels and fill:
	print("-f and -e are incompatible! Please decide if you want all sites or not")
	sys.exit(2)

# Get the header info
headerinfo = extractHeaderInfo(input)	#skips header, retains IDs in headerinfo

# Get the sample labels
IDs = []	#IDs holds all sample names without equal spacing (used in nexus)
IDs.append("reference ")
resultsequences = []	#resultsequences holds all complete sequences to be written to phylip or nexus file
resultsequences.append([])		#resultsequences[0] holds reference
for entry in headerinfo[1]:
	IDs.append(entry.replace("-",".")+" ")
	resultsequences.append([])
samplenames = fillUp(IDs)
	
linecounter = 0
print("\ngenerating phylip file with ",len(samplenames)-1," individuals")

# Go through the lines to get the genotypes
for line in input:
	site = line.strip('\n').split('\t')
	indel = False
	pos = int(site[1])
	
	# If missing positions should be filled up with Ns (-f specified, e.g. sites of low quality that were filtered out)
	if pos > (prev + 1) and fill:
		addLine=pos-(prev+1)
		individualcounter = 1
		for individual in site[9:]:
			resultsequences[individualcounter]+= "N" * addLine
			individualcounter += 1
		linecounter += addLine
		resultsequences[0] += "N" * addLine

	# site contains a deletion, replace by missing data
	if len(site[3])>1:
		indel=True

	else:
		alleles = site[4].split(",")  # if there are more than 1 alternative alleles, they are separated by commas
		for alt in alleles:
			if len(alt)>1 or '*' in alt:  # in case of an insertion
				indel=True
				break
	alternativeslist = []
	alternativeslist.append(site[3])
	if site[4] != ".":
		for entry in site[4].split(","):
			alternativeslist.append(entry)
	individualcounter = 1

	if not indel:
		for individual in site[9:]:
			gt = parse_gt(individual)
			if gt is None:
				resultsequences[individualcounter] += "N"
			else:
				resultsequences[individualcounter] += get_iupac(gt, alternativeslist)
			individualcounter += 1
		linecounter += 1
		resultsequences[0] += site[3]

	# If the site is an indel and noIndels is not specified, print as missing data (else not printed)
	elif not noIndels:
		for individual in site[9:]:
			resultsequences[individualcounter] += "N"
			individualcounter += 1
		linecounter += 1
		resultsequences[0] += site[3][:1]


	prev=int(site[1])

input.close()
	
writePhylipSequences(samplenames, resultsequences, output, writeref)
output.write("\n")
output.close()