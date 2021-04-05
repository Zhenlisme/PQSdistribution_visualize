##2020.07.24
"""
This script was writtern for calculating the G4 conservation value in all viruses.
"""

import argparse
import os,re,math,sys
from collections import defaultdict
import time
    
def G4detection(genomedict,Grun):
    complement={'A':"T","T":"A","G":"C","C":"G",
                "K":"M","M":"K","Y":"R","R":"Y","S":"S","W":"W",
                "B":"V","V":"B","H":"D","D":"H","N":"N","X":"X"}
    G4_dict = defaultdict(lambda:defaultdict(list))
    for ncnumber in genomedict:
        sequence=genomedict[ncnumber].replace('-','').upper()
        G4iter_positive=re.finditer('(?:G{%s}.{1,7}){3,}G{%s}'%(Grun,Grun), sequence)
        G4iter_negative=re.finditer('(?:C{%s}.{1,7}){3,}C{%s}'%(Grun,Grun), sequence)
        for g4 in G4iter_positive:
            start,end,g4seq=str(g4.start() + 1), str(g4.end()), g4.group()
            G4_dict[ncnumber]["+"].append([int(start),int(end),g4seq])
        for g4 in G4iter_negative:
            start,end,g4seq=str(g4.start() + 1), str(g4.end()), g4.group()
            g4seq=list(g4seq)
            g4seq.reverse()
            g4seq="".join([complement[i] for i in g4seq])
            G4_dict[ncnumber]["-"].append([int(start),int(end),g4seq])
    if G4_dict:
        return G4_dict
    else:
        print("No G%s-PQSs detected!"%Grun)
        sys.exit(0)

def intersect(location1, location2,slip=0):
    if location1[0]-slip > location2[1] or location1[1] < location2[0]-slip:
        return False
    else:
        return True
def locating(location,sequence,intbase,intgap):
    gapcount=intgap
    basecount=intbase
    for base in sequence:
        basecount+=1
        if base=="-":
            gapcount+=1
        if basecount-gapcount==location:
            start=basecount
            return (start,basecount,gapcount)

def location_correct(G4_dict,genome_dict):
    primary_correct=defaultdict(dict)
    for ncnumber in G4_dict.keys():
        extension_NC=defaultdict(list)
        genome_seq = genome_dict[ncnumber]
        genomelength=len(genome_seq.replace('-',''))
        for strand in ['+','-']:
            for NC_location in G4_dict[ncnumber][strand]:
                NCs,NCe=NC_location[0]-Args.elong,NC_location[1]+Args.elong
                NCs=NCs if NCs>=1 else 1
                NCe=NCe if NCe<=genomelength else genomelength
                extension_NC[NCs].append([NC_location[0]])
                extension_NC[NCe].append([NC_location[1]])
        extension_values=sorted(extension_NC.keys())
        gapcount,basecount=0,0
        for location in extension_values:
            corrected_location,basecount,gapcount=locating(location,genome_seq[basecount:],basecount,gapcount)
            for pqslist in extension_NC[location]:
                for pqssite in pqslist:
                    primary_correct[ncnumber][pqssite]=corrected_location
    return primary_correct
        
def G4SCI(aligned_genome,Grun):
    genome_dict={}
    with open(aligned_genome,'r') as F:
        for line in F:
            if line.startswith(">"):
                key=re.split('\s',line)[0].strip('>')
                genome_dict[key]=""
            else:
                genome_dict[key]+=line.rstrip()
                  
    G4_dict=G4detection(genome_dict,Grun)
    totalnumber = len(genome_dict.keys())
    NCL = location_correct(G4_dict,genome_dict)
    G4SCI_dict = defaultdict(lambda:defaultdict(dict))
    for ncnumber in G4_dict:
        for strand in ['+','-']:
            for NC_location in G4_dict[ncnumber][strand]:
                correct_location=(NCL[ncnumber][NC_location[0]],NCL[ncnumber][NC_location[1]])
                G4SCI_dict[ncnumber][strand][correct_location]=[1,NC_location]
                for strain_number in set(G4_dict.keys())-{ncnumber}:
                    for strain_location in G4_dict[strain_number][strand]:
                        strain_correct_location=(NCL[strain_number][strain_location[0]],NCL[strain_number][strain_location[1]])
                        if intersect(correct_location,strain_correct_location):
                            G4SCI_dict[ncnumber][strand][correct_location][0]+=1
                            break
    oplines=[]
    for ncnumber in G4SCI_dict.keys():
        for strand in G4SCI_dict[ncnumber]:
            for location in G4SCI_dict[ncnumber][strand]:
                g4sci = str(round(G4SCI_dict[ncnumber][strand][location][0] / totalnumber, 3))
                start, end, g4seq = G4SCI_dict[ncnumber][strand][location][1]
                oplines.append("\t".join([ncnumber,str(start),str(end),strand, str(totalnumber), str(g4sci), g4seq]) + "\n")
    return oplines

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="To get the location and conservation value of PQSs sequence. First-writtern@20190716 by Li Zhen.")
    parser.add_argument("-AS", "--aligned_sequence", type=str,help="The aligned genome sequence.")
    parser.add_argument("-e","--elong",default=0,type=int,help="Please type in the length of extension (nt).")
    parser.add_argument("-grun","--grun",type=int,help="Please input the length of grun.")
    parser.add_argument("-opf", "--opfile",type=str, help="The output file.")
    Args = parser.parse_args()
    with open(os.path.abspath(Args.opfile),'w') as F:
        oplines=G4SCI(os.path.abspath(Args.aligned_sequence),Args.grun)
        F.write("\t".join(["ID","start","end","strand","alignedcount","conservation","PQS"])+"\n")
        F.writelines(oplines)

