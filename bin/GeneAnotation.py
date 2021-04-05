from collections import defaultdict
import re,os
import argparse


class PQSinGene_coord():
    def __init__(self,gff,PQSfile=""):
        self.gff=gff
        self.PQSfile=PQSfile
    def anotate_utr(self):
        anotation_dir=defaultdict(lambda :defaultdict(list))
        with open(self.gff,'r') as F:
            for line in F:
                if line.startswith("##sequence-region"):
                    MatchKey,region_start,region_end=line.rstrip().split(" ")[1:]
                    anotation_dir[MatchKey]['region'] = [int(region_start),int(region_end)]
                    continue
                MatchKey=re.match(r'[A-Za-z]+?_?\d+\.?\d?',line,re.I)
                if MatchKey:
                    if line.split("\t")[2]=="CDS":
                        genename=re.findall(r'Parent=(.+?);',line,re.I)[0] if "Parent=" in line else \
                                 re.findall(r'Name=(.+?);',line,re.I)[0]
                        strandtype=line.split("\t")[6]
                        start,end=[int(i) for i in line.split("\t")[3:5]]
                        anotation_dir[MatchKey[0]][genename].append((start,end,'CDS',strandtype))
                        limit_left,limit_right=anotation_dir[MatchKey[0]]['region']
                        if strandtype == '+':
                            futr=(start-self.limit,start-1,'futr',strandtype) if start-self.limit>=limit_left \
                                  else (limit_left,start-1,'futr',strandtype)
                            tutr=(end+1,end+self.limit,'tutr',strandtype) if end+self.limit<=limit_right \
                                  else (end+1,limit_right,'tutr',strandtype)
                        else:
                            tutr=(start-self.limit,start-1,'tutr',strandtype) if start-self.limit>=limit_left \
                                  else (limit_left,start-1,'tutr',strandtype)
                            futr=(end+1,end+self.limit,'futr',strandtype) if end+self.limit<=limit_right \
                                  else (end+1, limit_right,'futr',strandtype)
                        anotation_dir[MatchKey[0]][genename].extend([futr,tutr])
        return anotation_dir
    def anotate_noncoding(self):
        anotation_dir = defaultdict(lambda: defaultdict(list))
        opanodir=defaultdict(lambda :defaultdict(list))
        with open(self.gff,'r') as F:
            for line in F:
                if line.startswith("##sequence-region"):
                    MatchKey,region_start,region_end=line.rstrip().split(" ")[1:]
                    anotation_dir[MatchKey]['region'] = [int(region_start),int(region_end)]
                    continue
                MatchKey=re.match(r'[A-Za-z]+?_?\d+\.?\d?',line,re.I)
                if MatchKey:
                    if line.split("\t")[2]=="CDS":
                        genename=re.findall(r'Parent=(.+?);',line,re.I)[0] if "Parent=" in line else \
                                 re.findall(r'Name=(.+?);',line,re.I)[0]
                        strandtype=line.split("\t")[6]
                        start,end=[int(i) for i in line.split("\t")[3:5]]
                        anotation_dir[MatchKey[0]][genename].append((start,end,'CDS',strandtype))
        for ncnumber in anotation_dir:
            if not [i for i in anotation_dir[ncnumber] if i!="region"]:
                continue
            left_limit, right_limit = anotation_dir[ncnumber]['region']
            positive_location=[[location[0:2] for location in anotation_dir[ncnumber][gene] if location[-1]=="+"]
                                              for gene in anotation_dir[ncnumber] if gene!="region"]
            negative_location=[[location[0:2] for location in anotation_dir[ncnumber][gene] if location[-1]=="-"]
                                              for gene in anotation_dir[ncnumber] if gene!="region"]
            positive_location = sorted([b for a in positive_location for b in a],key=lambda X:(X[0],X[1]))
            negative_location = sorted([b for a in negative_location for b in a],key=lambda X:(X[0],X[1]))
            positive_noncoding = [(positive_location[i][1] + 1, positive_location[i + 1][0] - 1)
                                 for i in range(len(positive_location)) if len(positive_location) > i + 1 and
                                 positive_location[i + 1][0] > positive_location[i][1]+2]
            if positive_location:
                pemax=sorted([i[1] for i in positive_location])[-1]
                if positive_location[0][0]>left_limit:
                    positive_noncoding.insert(0,(left_limit,positive_location[0][0]-1))
                if pemax<right_limit:
                    positive_noncoding.append((pemax+1,right_limit))
            else:
                positive_noncoding.insert(0,(left_limit,right_limit))
            negative_noncoding= [(negative_location[i][1] + 1, negative_location[i + 1][0] - 1)
                                 for i in range(len(negative_location)) if len(negative_location) > i + 1 and
                                 negative_location[i + 1][0] > negative_location[i][1]+2]
            if negative_location:
                nemax=sorted([i[1] for i in negative_location])[-1]
                if negative_location[0][0]>left_limit:
                    negative_noncoding.insert(0,(left_limit,negative_location[0][0]-1))
                if nemax<right_limit:
                    negative_noncoding.append((nemax+1,right_limit))
            else:
                negative_noncoding.insert(0,(left_limit,right_limit))
            positive_coding = [(positive_noncoding[loc][1]+1,positive_noncoding[loc+1][0]-1) for loc in range(len(positive_noncoding)-1)]
            if positive_noncoding:
                if positive_noncoding[0][0]>left_limit:
                    positive_coding.insert(0,(left_limit,positive_noncoding[0][0]-1))
                if positive_noncoding[-1][-1]<right_limit:
                    positive_coding.append((positive_noncoding[-1][1]+1,right_limit))
            else:
                positive_coding.append((left_limit,right_limit))
            negative_coding = [(negative_noncoding[loc][1]+1,negative_noncoding[loc+1][0]-1)
                               for loc in range(len(negative_noncoding)-1)]
            if negative_noncoding:
                if negative_noncoding[0][0]>left_limit:
                    negative_coding.insert(0,(left_limit,negative_noncoding[0][0]-1))
                if negative_noncoding[-1][-1]<right_limit:
                    negative_coding.append((negative_noncoding[-1][1]+1,right_limit))
            else:
                negative_coding.append((left_limit,right_limit))
            positive_noncoding=[(loc[0],loc[1],"_".join(['noncoding',str(positive_noncoding.index(loc))]),'+')
                                for loc in positive_noncoding]
            negative_noncoding=[(loc[0],loc[1],"_".join(['noncoding',str(negative_noncoding.index(loc))]),'-')
                                for loc in negative_noncoding]
            positive_coding=[(loc[0],loc[1],"_".join(['coding',str(positive_coding.index(loc))]),'+')
                             for loc in positive_coding]
            negative_coding=[(loc[0],loc[1],"_".join(['coding',str(negative_coding.index(loc))]),'-')
                             for loc in negative_coding]
            opanodir[ncnumber]['non_coding'].extend(positive_noncoding)
            opanodir[ncnumber]['non_coding'].extend(negative_noncoding)
            opanodir[ncnumber]['coding'].extend(positive_coding)
            opanodir[ncnumber]['coding'].extend(negative_coding)
        return opanodir
    def anotate_nongene(self):
        anotation_dir = defaultdict(lambda: defaultdict(list))
        with open(self.gff, 'r') as F:
            for line in F:
                if line.startswith("##sequence-region"):
                    MatchKey, region_start, region_end = line.rstrip().split(" ")[1:]
                    anotation_dir[MatchKey]['region'] = [int(region_start), int(region_end)]
                    continue
                MatchKey = re.match(r'[A-Za-z]+?_?\d+\.?\d?', line, re.I)
                if MatchKey:
                    if line.split("\t")[2] == "CDS":
                        genename = re.findall(r'Parent=(.+?);', line, re.I)[0] if "Parent=" in line else \
                            re.findall(r'Name=(.+?);', line, re.I)[0]
                        genename = genename.replace("gene-","")
                        strandtype = line.split("\t")[6]
                        start, end = [int(i) for i in line.split("\t")[3:5]]
                        anotation_dir[MatchKey[0]][genename].append((start, end, 'CDS', strandtype))
        for ncnumber in anotation_dir:
            if not [i for i in anotation_dir[ncnumber] if i!="region"]:
                continue
            left_limit, right_limit = anotation_dir[ncnumber]['region']
            codings=[[location[0:2] for location in anotation_dir[ncnumber][gene]]
                     for gene in anotation_dir[ncnumber] if gene != "region"]
            codings=sorted([b for a in codings for b in a], key=lambda X: (X[0], X[1]))
            # noncodings = [(codings[i][1] + 1, codings[i + 1][0] - 1) for i in range(len(codings))
            #              if len(codings) > i + 1 and codings[i + 1][0] > codings[i][1]]
            noncodings=[]
            for i in range(len(codings)):
                if len(codings) > i + 1 and codings[i + 1][0] > codings[i][1]:
                    maxvalue=max([a[1] for a in codings[:i+1]])
                    if maxvalue<codings[i + 1][0]:
                        noncodings.append([maxvalue+1, codings[i + 1][0] - 1])
            if codings:
                pemax=max([i[1] for i in codings])
                if codings[0][0]>left_limit:
                    noncodings.insert(0,(left_limit,codings[0][0]-1))
                if pemax<right_limit:
                    noncodings.append((pemax+1,right_limit))
            else:
                noncodings.insert(0,(left_limit,right_limit))
            anotation_dir[ncnumber]['noncoding']=[i for i in noncodings if i[0]<i[1]]
        return {nc:[(b[0],b[1],a,b[-1]) if a!="noncoding" else (b[0],b[1],a) for a in anotation_dir[nc]
                    for b in anotation_dir[nc][a] if a!="region"] for nc in anotation_dir }
    def creatanotation(self,optfile,limit=150):
        self.opf = optfile
        self.limit = limit
        if self.limit=="noncoding":
            anotation=self.anotate_noncoding()
        elif isinstance(self.limit,int):
            anotation=self.anotate_utr()
        else:
            print("wrong parameter!")
            exit(0)
        opline=[]
        for nc in anotation:
            for gene in anotation[nc]:
                if gene=="region":
                    continue
                else:
                    [opline.append([nc,str(location[0]),str(location[1]),gene,location[2],location[3]])
                     for location in anotation[nc][gene]]
        opline=sorted(opline,key=lambda X:(X[0],X[-1],int(X[1]),X[3]))
        with open(self.opf, 'w') as F:
            F.writelines(["\t".join(i)+"\n" for i in opline])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="To get the location and conservation value of PQSs sequence. First-writtern@20210405 by Li Zhen.")
    parser.add_argument("-a", "--anotation", default="",type=str, help="Please input the anotation file")
    parser.add_argument("-opf", "--opfile",type=str, help="The output file.")
    Args = parser.parse_args()
    nongene=PQSinGene_coord(gff=os.path.abspath(Args.anotation)).anotate_nongene()
    oplines=[]
    for nc in nongene:
        genelist=sorted(nongene[nc],key=lambda x:x[0])
        i = 1
        for gene in genelist:
            if len(gene)==3:
                ano="".join(['noncoding','(',str(i),')'])
                i+=1
            else:
                ano="".join([gene[2],'(',gene[-1],')'])
            oplines.append("\t".join([nc, str(gene[0]), str(gene[1]), ano])+"\n")
    with open(os.path.abspath(Args.opfile),'w') as F:
        F.writelines(oplines)

