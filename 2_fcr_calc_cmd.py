#!/usr/bin/python
import sys
import math
# Load as a map the csv with PMID in first column and its matched Journal Citation Rate in the second
jcrmap = dict()
with open("./intermed_data/jcrmap.csv") as file:
    for line in file:
        info = line.split(",")
        pmid = info[0]
        jcr = float(info[1])
        jcrmap[pmid] = jcr

# Load citation tree tsv. PMID in first column cites the PMIDs in the second column (cited articles are pipes-separated)
citingmap = dict()
i = 0
with open("./intermed_data/mergeCites.tsv") as file:
    for line in file:
        info = line.split("\t")
        i += 1
        citing = info[0].strip(" \"\r\n")
        if citing == "pmid":
            continue
        #if i % 1000000 == 0:
        #    print(i)
        cited = info[1].strip(" \"\r\n")
        citedList = cited.split("|")
        citingmap[citing] = [x.strip(" \"\r\n") for x in citedList]
    del i

# Stdin feeds cited-by article tree (tsv).  PMID in first column cited by PMIDs in second column (citing articles are pipes separated)
i = 0
for line in sys.stdin:
    line = line.strip(" \"\r\n")
    if line == "":
        continue
    i +=1
    pmids = line.split("\t")
    # Cited article is our article of interest
    cited = pmids[0].strip(" \"\r\n")
    # Find all the articles citing it as the first step for building co-citation network
    citing = pmids[1].strip(" \"\r\n").split("|")
    # Find all articles cited by the articles citing the article of interest. This is the co-citation network.
    network = list()
    for pmid in citing:
        try:
            network.extend(citingmap[pmid])
        except:
            pass
    if len(network) == 0:
        continue
    jcrs = list()
    # De-duplicate co-citation network
    network = list(set(network))
    # Find JCRs for all articles in co-citation network
    for x in network:
        try:
            jcrs.append(jcrmap[x.strip(" \"\r\n")])
        except:
            pass
    if len(jcrs) == 0:
        continue
    # Field Citaiton Rate is the average of all Journal Citation Rates for the co-citation network articles.
    fcr = math.fsum(jcrs) / len(jcrs)
    print(cited + "\t" + str(fcr))

