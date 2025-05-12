#! /usr/bin/env python

import os
import sys
import re
from collections import defaultdict

es  = defaultdict(list) #retrieve exon start give an exon id
ee  = defaultdict(list) #retrieve exon end given an exon id
eg  = defaultdict(list) #retrieve exon ids given an gene id
ese = defaultdict(list)	#retrieve sequences given an exon id
ec  = defaultdict(list) #retrieve chromosome information for an exon id
est = defaultdict(list) #retrieve strand Information for an exon id
en  = defaultdict(list) #retrieve exon number for an exon id
et  = defaultdict(list)

class ParseGTF():
    def __init__(self, gtf):
        f=open(gtf)
        line = f.readline()
        while line:
            if not (re.match(r'#', line)):
                exn=re.findall(r'exon_id "(.*?)";',line)
                if(len(exn)>=1):
                    es[exn[0]].append(line.split()[3])
                    ee[exn[0]].append(line.split()[4])
                    ec[exn[0]].append(line.split()[0])
                    est[exn[0]].append(line.split()[6])		
                    tr=re.findall(r'transcript_id "(.*?)";',line)
                    gn = re.findall(r'gene_id "(.*?)";',line)
                    enu = re.findall(r'exon_number "(.*?)";',line)
                    et[exn[0]].append(tr[0])
                    eg[gn[0]].append(exn[0])
                    eg[tr[0]].append(exn[0])
                    en[exn[0]].append(enu[0])
            line=f.readline()
		
    def getExons(self, key):
        return eg[key]
	
    def getStart(self, key):
        return int(es[key][0])

    def getEnd(self, key):
        return int(ee[key][0])

    def getChromosome(self, key):
        return ec[key]

    def getStrand(self, key):
        return est[key]

    def getExonNumber(self, key):
        return int(en[key][0])

    def getTranscript(self, key):
        return et[key][0]

    def getAllGenes(self):
        i=0
        genes = []
        for key in eg.keys():
            if not (key.find("ENSG")):
                genes.append(key)
                i=i+1
        return genes

    def getTranscriptPosition(self, key):
        allexons = eg[key]
        startpoint=0
        endpoint=0
        if(len(allexons)>0):
            startpoint=int(ee[allexons[0]][0])
            endpoint=0
            for e in allexons:
                en=int(ee[e][0])
                if(en>endpoint):
                    endpoint=en
        return startpoint,endpoint

    def getAllTranscripts(self):
        i=0
        transcripts = []
        for key in eg.keys():
            if not (key.find("ENST")):
                transcripts.append(key)
                i=i+1
        return transcripts
