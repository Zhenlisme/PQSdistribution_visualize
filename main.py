import argparse,os

parser = argparse.ArgumentParser(description="To get the location and conservation value of PQSs sequence. First-writtern@20210406 by Li Zhen.")
parser.add_argument("-AS", "--aligned_sequence", type=str,help="The aligned genome sequence.")
parser.add_argument("-a", "--anotation", default="",type=str, help="Please input the anotation file")
parser.add_argument("-e","--elong",default=0,type=int,help="Please type in the length of extension (nt).")
parser.add_argument("-grun","--grun",type=int,help="Please input the length of grun.")
parser.add_argument("-text","--text", type=str,help="Anotation to the circle plot.")
parser.add_argument("-opf","--opfile", type=str,help="output file name.")
Args = parser.parse_args()

##Users uploaded
Alignfile=os.path.abspath(Args.aligned_sequence)
gffile=os.path.abspath(Args.anotation)

####To be deleted
anobed=os.path.abspath('.Ano.bed')
G4sci=os.path.abspath('.g4sci')
G4_coord=os.path.abspath('.g4coord')

###Results
circlplot=os.path.abspath(Args.opfile)+'.jpg'

##AnoBed
os.system('python3 %s -a %s -opf %s'%(os.path.abspath("bin/GeneAnotation.py"),gffile,anobed))

##G4SCI
os.system('python3 %s -AS %s -e %d -grun %d -opf %s'%(os.path.abspath("bin/conservation.py"),Alignfile,Args.elong,Args.grun,G4sci))

##G4_coord
os.system('python3 %s -a %s -scif %s -opf %s'%(os.path.abspath("bin/G4_coord.py"),anobed,G4sci,G4_coord))

##draw plot
PQStype="".join(["G",str(Args.grun),"-PQS"])
os.system('Rscript %s %s %s %s %s'%(os.path.abspath("bin/circlplot_simplified.r"),G4_coord,PQStype,Args.text,circlplot))

#delete files
os.system("rm -f %s %s %s"%(anobed,G4sci,G4_coord))