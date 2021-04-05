##2020.07.24
"""
This script was writtern for calculating the GNG SCI in all viruses.
"""

import argparse
import os,re,math
from collections import defaultdict
import time

def PQScoord(G4SCIfile,anotation_bed,opfile):
    G4SCI_dict=defaultdict(list)
    anotation_pqs=[]
    with open(G4SCIfile,'r') as F:
        for line in F:
            ncnumber,start,end,strand,strain_count,g4sci,g4seq=line.rstrip().split("\t")
            G4SCI_dict[ncnumber].append((start,end,g4sci,strain_count,strand))
    # already_list=[]
    with open(anotation_bed,'r') as F:
        for line in F:
            ncnumber, start, end,anotat = line.rstrip().split("\t")
            coord_list=[]
            coord_list.append(start)
            g4sci_pqscount=[]
            if ncnumber not in G4SCI_dict:
                continue
            for pqs in G4SCI_dict[ncnumber]:
                strain_count=pqs[3]
                # already_list.append(ncnumber)
                if int(pqs[0])<=int(end) and int(pqs[0])>=int(start):
                    coord_list.extend([str(int(pqs[0])-1),str(int(pqs[0])),str(int(pqs[0])+10),str(int(pqs[0])+11)])
                    pqscount = '1' if pqs[-1] == "+" else '-1'
                    g4sci_pqscount.append([pqs[2],pqscount])
            if int(coord_list[-1])<=int(end):
                coord_list.append(end)
            else:
                coord_list.append(coord_list[-1])
            anotation_pqs.extend(["\t".join((ncnumber,coord_list[i*2],coord_list[i*2+1],'0','0',anotat,strain_count))+"\n"
                    if i%2==0 else "\t".join((ncnumber,coord_list[i*2],coord_list[i*2+1],g4sci_pqscount[int((i-1)/2)][0],
                    g4sci_pqscount[int((i-1)/2)][1],anotat,strain_count))+"\n" for i in range(int(len(coord_list)/2))])
    with open(opfile,'w') as F:
        anotation_pqs=[i for i in anotation_pqs if int(i.split("\t")[2])>int(i.split("\t")[1])]
        F.writelines(anotation_pqs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="To get the location and anotation of a sequence. First-writtern@20190716 by ZhenLi.")   
    parser.add_argument("-a", "--anotation", default="",type=str, help="Please input the anotation file")
    parser.add_argument("-opf", "--opfile",type=str, help="The output file.")
    parser.add_argument("-scif","--g4sf",default="",type=str,help="Please type in the g4sci file.")
    Args = parser.parse_args()
    PQScoord(os.path.abspath(Args.g4sf),os.path.abspath(Args.anotation),os.path.abspath(Args.opfile))
