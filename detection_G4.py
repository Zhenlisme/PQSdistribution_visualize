import re,argparse,os

def G4finder(genomefile,opfile,Grun):
    genome_dict={}
    with open(genomefile, 'r') as F:
        for line in F:
            if line.startswith(">"):
                key = line.split(" ")[0].strip(">\n")
                genome_dict[key] = ""
            else:
                genome_dict[key] += line.rstrip().upper()
    oplines=[]
    for ncnumber in genome_dict:
        G4iter_positive=re.finditer('(?:G{%s}.{1,7}){3,}G{%s}'%(Grun,Grun), genome_dict[ncnumber])
        G4iter_negative=re.finditer('(?:C{%s}.{1,7}){3,}C{%s}'%(Grun,Grun), genome_dict[ncnumber])
        for g4 in G4iter_positive:
            start,end,g4seq=str(g4.start() + 1), str(g4.end()), g4.group()
            oplines.append(",".join([ncnumber,start,end,'G'*Grun,g4seq])+"\n")
        for g4 in G4iter_negative:
            start,end,g4seq=str(g4.start() + 1), str(g4.end()), g4.group()
            oplines.append(",".join([ncnumber,start,end,'C'*Grun,g4seq])+"\n")
    with open(opfile,'w') as F:
        F.writelines(oplines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="To get the location of PQS. First-writtern@20190716 by ZhenLi.")
    parser.add_argument("-f", "--genomefile", type=str,required=True,help="Please type in the file of genome file.")
    parser.add_argument("-g","--g4type",type=int,choices=[2,3],required=True,help="Please type in the PQS type.")
    parser.add_argument("-o", "--opfile", required=True, type=str, help="The output file name.")
    Args = parser.parse_args()
    G4finder(os.path.abspath(Args.genomefile),os.path.abspath(Args.opfile),int(Args.g4type))

