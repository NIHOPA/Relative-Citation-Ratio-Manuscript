#! python
# -*- coding: utf-8 -*-

"""
Created in 2014

@author: B. Ian Hutchins
"""
import sys
import math
import json
import collections
methods = ["sibling", "parent", "child"]

# Cast variables to floats with error handling
def floatNaN(f):
    try:
        return(float(f))
    except (ValueError, TypeError):
        return(float('nan'))

# Cast variables to strings with error handling
def strEmpty(s):
    try:
        if s is None:
            return('')
        else:
            return(str(s))
    except:
        return('')

# Function to load precleaned data from JSON files
def loadCiteDataCleaned(years):
    numColnames = ["ISI_ItemNumberID", "ISI_ItemID", "ISI_PubYear"]
    charColnames = ["ISI_References"]
    charColnames.extend(["X" + str(x) for x in years])
    charColnames.extend(["ISI_FullSourceTitle"])
    
    data = dict()
    
    for i in years:
        print("Getting year " + str(i))
        dataYr = json.load(open((settings["trDataPath"] + str(i) + ".json")))
        print("Extract and append " + str(i))
        
        for j in range(len(dataYr["ISI_ItemID"])):
            trid = dataYr["ISI_ItemID"][j]
            data[trid] = dict()
            
            for k in numColnames:
                data[trid][k] = floatNaN(dataYr[k][j])
            
            for k in charColnames:
                data[trid][k] = strEmpty(dataYr[k][j])
    
    return(data)

# Function to construct an article ID lookup table
def makeT9tridTable(data):
    table = collections.defaultdict(list)
    for i in list(data.keys()):
        if data[i]["ISI_ItemNumberID"] != float('nan'):
            table[data[i]["ISI_ItemNumberID"]].append(i)
    # Delete any erroneous records that have '0' as an article ID
    del(table[0])
    return(table)

# Function to load a table of journal citation rates
def loadJournalJifs(yrs):
    jjdict = dict()
    for i in yrs:
        dataYr = json.load(open((settings["jifsDataPath"] + str(i) + ".json")))
        jjdict[i] = dict()
        for j in range(len(dataYr["SRC_TITLE"])):
            jjdict[i][dataYr["SRC_TITLE"][j]] = floatNaN(dataYr["JCR_IF"][j])
    return(jjdict)

# Function to lookup the IDs of articles cited by a specified article of interest
def getChildT9s(trid, data):
    t9s = data[trid]["ISI_References"].split("|")
    if t9s.count("") > 0:
        for i in range(t9s.count("")):
            t9s.remove("")
    
    if len(t9s) == 0:
        return(t9s)
    
    t9s = [floatNaN(x) for x in t9s]
    if t9s.count(float('nan')) > 0:
        for i in range(t9s.count(float('nan'))):
            t9s.remove(float('nan'))
    
    return(t9s)

# Function to lookup the IDs of articles citing a specified article of interest
def getParentT9s(trid, data):
    t9Colnames = (["X" + str(x) for x in list(range(settings["startYear"],settings["endYear"]))])
    t9s = list()
    for c in t9Colnames:
        temp = data[trid][c].split("|")
        if temp.count("") > 0:
            for i in range(temp.count("")):
                temp.remove("")
        t9s.extend(temp)
        
    t9s = [floatNaN(x) for x in t9s]
    
    if len(t9s) == 0:
        return(list())
    elif t9s.count(float('nan')) > 0:
        for i in range(t9s.count(float('nan'))):
            t9s.remove(float('nan'))
    return(t9s)

# Function to lookup the IDs of articles co-cited with a specified article of interest
def getSiblingT9s(trid, data, table):
    # Get IDs of citing articles
    t9sParent = getParentT9s(trid, data)
    if len(t9sParent) > 0:
        tridsParent = list()
        for i in t9sParent:
            tridsParent.extend(table[i])
        t9sSiblings = list()
        # Look up other articles cited along with the article of interest
        for x in tridsParent:
            t9sSiblings.extend(getChildT9s(x, data))
        t9sSiblings = [floatNaN(x) for x in t9sSiblings]
        if t9sSiblings.count(float('nan')) > 0:
            for i in range(t9sSiblings.count(float('nan'))):
                t9sSiblings.remove(float('nan'))
        
        t9sSiblings = list(set(t9sSiblings))
        return(t9sSiblings)
    else:
        return(list())

# Function to find the articles in the specified level of the citation network
# (citing, cited, or co-cited), look up the journal citation rates, and average
# these to generate the Field Citation Rate, for the specified article of interest.
# If ret = 'default', then the Field Citation Rate will be returned, otherwise
# the journal citation rates will be jointed with '|' and returned as a string
def cocitationJifs(trid, journJifs, data, t9TRIDtable, method="sibling", ret="default"):
    # Generate list of article IDs for the specified level of the citation network
    # 'parent' = citing articles, 'child' = cited articles, and 'sibling' = co-cited articles
    t9s = list()
    if method == "parent":
        t9s = getParentT9s(trid, data)
    elif method == "sibling":
        t9s = getSiblingT9s(trid, data, t9TRIDtable)
    else:
        t9s = getChildT9s(trid, data)
    for i in range(t9s.count("")):
        t9s.remove("")
    
    # Clean 0s from the list
    t9s = list(set(t9s))
    if t9s.count(0) > 0:
        t9s.remove(0)
    
    # Translate T9 article ID numbers to accession numbers
    trids = list()
    for i in t9s:
        trids.extend(t9TRIDtable[i])
    
    if trid not in trids:
        trids.append(trid)
    
    # Look up journal and publication year for each article in the network
    # and find the corresponding journal citation rate
    jInfo = list()
    yInfo = list()
    jInfo = [data[x]["ISI_FullSourceTitle"] for x in trids]
    yInfo = [data[x]["ISI_PubYear"] for x in trids]
    jifs = list()
    for n in list(range(len(jInfo))):
        try:
            jifs.append(journJifs[int(yInfo[n])][jInfo[n]])
        except:
            pass
    
    # 
    jifs = [floatNaN(x) for x in jifs]
    for i in range(jifs.count(float('nan'))):
        jifs.remove(float('nan'))
    
    if len(jifs) > 0:
        if ret == "default":
            return(math.fsum(jifs) / len(jifs))
        else:
            return("|".join([str(x) for x in jifs]))
    else:
        return(float('nan'))

# Give the path to the configuration file as an argument to the script
if __name__=='__main__':
    configPath = sys.argv[1]
    
    settings = dict()
    
    settings["trDataPath"] = "./"
    settings["jifsDataPath"] = "./"
    settings["matchedDataPath"] = "./"
    settings["outputPath"] = "./"
    
    with open(configPath) as file:
        for line in file:
            info = line.split("=")
            settings[info[0]] = info[1].strip("\"\r\n")
    
    try:
        settings["startYear"] = int(settings["startYear"])
    except:
        settings["startYear"] = 2002
    
    try:
        settings["endYear"] = int(settings["endYear"]) + 1
    except:
        settings["endYear"] = 2013
    
    # Load table of journal citation rates, 'jj'
    jj = loadJournalJifs(list(range(settings["startYear"], settings["endYear"])))
    # Load citaiton data, 'data'
    data = loadCiteDataCleaned(list(range(settings["startYear"], settings["endYear"])))
    # Load table to convert T9 article IDs to Web of Science accession numbers
    t9table = makeT9tridTable(data)
    # Data files are split by year
    espaYears = list(range(settings["startYear"], settings["endYear"] - 1))
    # For this script, calculate Field Citation Rates for all three network levels
    for m in methods:
        for e in espaYears:
            # Read article IDs matched to Web of Science accession numbers 
            filetext = settings["matchedDataPath"] + str(e) + ".txt"
            matchFile = open(filetext, "r")
            match = list()
            for line in matchFile:
                match.append(line)
            
            matchFile.close()
            match = [x.strip() for x in match]
            match = [floatNaN(x) for x in match]
            
            # Open output file and write header
            outPath = settings["outputPath"] + m + " " + str(e) + ".txt"
            outfile = open(outPath, "a")
            outfile.write("TRID,FCR\n")
            # Calculate Field Citaiton Rates and write to file
            for i in range(len(match)):
                try:
                    outstring = str(int(match[i])) + "," + str(cocitationJifs(match[i], jj, data, t9table, method=m, ret="default")) + "\n"
                    outfile.write(outstring)
                    print(str(int(match[i])))
                except:
                    print("Exception")
                    outstring = str(int(match[i])) + "," + "None" + "\n"
                    outfile.write(outstring)
            
            
            outfile.close()
