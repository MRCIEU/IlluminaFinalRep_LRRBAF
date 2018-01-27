#!/usr/bin/python

import sys
import math
import numpy
import multiprocessing

# Script to estimate LRR and BAF using Final Report
# Data from external provider in Illumina Final Report with column headings:
#        SNP Name,Sample ID,Allele1 - Top,Allele2 - Top,GC Score,X,Y,X Raw,Y Raw

def getcentroids(data):
    return [numpy.median(data[0]),numpy.median(data[1])]

def calculateLRR(R,theta,centroids):
    try:
        if theta > centroids[1][1]:
            LRR = math.log(R/(centroids[1][0]+((centroids[2][0]-centroids[1][0])*abs((theta-centroids[1][1])/(centroids[1][1]-centroids[2][1])))),2)
        elif theta < centroids[1][1]:
            LRR = math.log(R/(centroids[1][0]+((centroids[0][0]-centroids[1][0])*abs((centroids[1][1]-theta)/(centroids[1][1]-centroids[0][1])))),2)  # should this be log2??
        else:
            LRR = math.log(R/centroids[1][0],2)
        return LRR
    except:
        return "NA"

def calculateBAF(R,theta,centroids):
    try:
        if theta > centroids[1][1]:
            BAF = 0.5+0.5*((theta-centroids[1][1])/(centroids[2][1]-centroids[1][1]))
        elif theta < centroids[1][1]:
            BAF = 0.5-0.5*((centroids[1][1]-theta)/(centroids[1][1]-centroids[0][1]))
        else:
            BAF = 0.5
        if BAF>1.0: BAF=1.0
        if BAF<0.0: BAF=0.0
        return BAF
    except:
        return "NA"

def calculatedata(data):
    genotypes = ["","",""]
    bygenotype = {}
    centroids = []
    result = ""
    p = []
    for line in data:
        p.append(line[2]+line[3])
    if not len(set(p)) == 3:
        return data[0][0], "# of Genotypes != 3, i.e. 2 Hom + 1 Het"
    try:
        for line in data:
            if line[2]+line[3] in genotypes:
                try:
                    bygenotype[line[2]+line[3]][0].append(float(line[5])+float(line[6]))
                    bygenotype[line[2]+line[3]][1].append(2.0/math.pi*math.atan(float(line[6])/float(line[5])))
                except:
                    null = 0
            else:
                if line[2] in ["A","C","G","T"] and line[3] in ["A","C","G","T"]:
                    if line[2] != line[3]:
                        genotypes[1] = line[2]+line[3]
                    elif genotypes[0] == "":
                        genotypes[0] = line[2]+line[3]
                    else:
                        genotypes[2] = line[2]+line[3]
                    bygenotype[line[2]+line[3]] = [0.0,0.0]
                    try:
                        bygenotype[line[2]+line[3]][0] = [float(line[5])+float(line[6])]
                        bygenotype[line[2]+line[3]][1] = [2.0/math.pi*math.atan(float(line[6])/float(line[5]))]
                    except:
                        null = 0
        for genotype in genotypes:
            centroids.append(getcentroids(bygenotype[genotype]))
        if centroids[2][1]<centroids[0][1]:
            tmp = centroids[2]
            centroids[2] = centroids[0]
            centroids[0] = tmp
        for line in data:
            if float(line[5])>0:
                LRR = calculateLRR(float(line[5])+float(line[6]),2.0/math.pi*math.atan(float(line[6])/float(line[5])),centroids)
            else:
                LRR = "NA"
            if float(line[5])>0:
                BAF = calculateBAF(float(line[5])+float(line[6]),2.0/math.pi*math.atan(float(line[6])/float(line[5])),centroids)
            else:
                BAF = "NA"
            line.append(str(LRR))
            line.append(str(BAF))
            result += "\t".join(line)+"\n"
        return result, None
    except Exception, e:
        return data[0][0], str(e)


if __name__ == '__main__':
    with open("./LRR/ALSPAC.raw", "r") as fin, \
         open("./LRR/ALSPAC.LRR", "w") as fo, \
         open("./LRR/LRR.out", "w") as out:
        X, Y  = 5, 6 # X/Y position in Final Report
        idx   = 0
        Batch = {}
        MP    = multiprocessing.Pool()
        for record in fin:
            record = record.strip().split("\t")
            if not Batch.has_key(record[0]):
                if len(Batch) == 2500:
                    idx += 1
                    print "Processing: Batch", idx

                    Data = MP.map(calculatedata, [ Batch[x] for x in sorted(Batch) ] )
                    for p in Data:
                        if p[1]:
                            out.write("\t".join([ p[0], p[1] ]) + "\n"); continue
                        fo.write(p[0])
                    Batch = {}
                Batch[record[0]] = []
            Batch[record[0]].append( record )
