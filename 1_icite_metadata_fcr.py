import collections
import csv
import json
import gc
import threading
import math


# Comma separated with quotes
# Headers = pmid, doi, year, date, journal, issn, cites, res, jcr
metaPath1 = "./start_data/medline16_meta_parsed.tsv"
metaPath2 = "./start_data/notMedline15_meta_parsed_noNA.csv"
# Tab separated with header = citing, cited
solrCitePath = "./start_data/solr_citations.tsv"

citesOutPath = "./intermed_data/mergeCites"
citedByOutPath = "./intermed_data/mergeCitedBy"
skippedOutPath = "./intermed_data/mergeSkipped.json"
metaOutPath = "./intermed_data/mergeMeta"


def floatNaN(f):
    try:
        return(float(f))
    except (ValueError, TypeError):
        return(float('nan'))

def strEmpty(s):
    try:
        if s is None:
            return('')
        else:
            return(str(s))
    except:
        return('')

def cleanString(s):
    try:
        if s is None:
            return('')
        # Get rid of zero-width spaces
        else:
            s = str(s)
            s2 = ''.join(['' if x == '\u200b' else x for x in s])
            s = ''.join([c if ord(c) > 31 and ord(c) != 127 else ' ' for c in s2])
            #s2 = ''.join([' ' if x == '\u00A0' else x for x in s])
            return(s)
    except:
        return('')


cites = collections.defaultdict(list)
citedBy = collections.defaultdict(list)
incJCR = list()
meta = collections.defaultdict(dict)
skipped = list()

# Split and clean row text, add data to data stores
def processRow(row):
    rowClean = [cleanString(x.strip("\"\'\r\n\t ")) if isinstance(x, str) else x for x in row]
    if len(rowClean) != 9 and len(rowClean) != 16:
        skipped.append(rowClean[0])
    if int(rowClean[8]) == 1:
        incJCR.append(int(rowClean[0]))
    meta[int(rowClean[0])] = {"doi":rowClean[1], "year":int(rowClean[2]), "journal":rowClean[4], "issn":rowClean[5], "res":int(rowClean[7]), "jcr":int(rowClean[8])}

print("Loading MedLine metadata")
# Headers = pmid, doi, year, date, journal, issn, cites, res, jcr
i = 0
with open(metaPath1, 'rb') as file:
    reader = csv.reader(file, quoting=csv.QUOTE_NONE, delimiter="\t")
    for row in reader:
        i += 1
        if i == 1:
            continue
        processRow(row)
        if i % 1000000 == 0:
            print(str(i))
    del i

print("Loading PubMed-Not-MedLine metadata")
i = 0
with open(metaPath2, 'rb') as file:
    reader = csv.reader(file, quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        i += 1
        if i == 1:
            continue
        processRow(row)
        if i >= 10:
            pass
    del i



print("Loading citation data")
i = 0
with open(solrCitePath, 'r') as file:
    for line in file:
        i += 1
        # Skip header
        if i == 1:
            continue
        lineData = line.split("\t")
        citingPmid = int(float(cleanString(lineData[0].strip("\"\'\r\n\t "))))
        citedPmid = int(float(cleanString(lineData[1].strip("\"\'\r\n\t "))))
        cites[citingPmid].append(citedPmid)
        citedBy[citedPmid].append(citingPmid)
        if i >= 1000:
            pass
        if i % 3000000 == 0:
            print("Loading citation data " + str(i))
    del i


print("Deduplicating citations")
# Deduplicate citation lists
for k in list(cites.keys()):
    temp = list(set(cites[k]))
    if k in temp:
        temp.remove(k)
    cites[k] = temp

for k in list(citedBy.keys()):
    temp = list(set(citedBy[k]))
    if k in temp:
        temp.remove(k)
    citedBy[k] = temp


# Calculate citation counts
print("Calculating citation counts")
for pmid in list(citedBy.keys()):
    meta[pmid]['citation_count'] = len(citedBy[pmid])

# Calculate number of references per paper
print("Calculating ref counts")
for pmid in list(cites.keys()):
    meta[pmid]['ref_count'] = len(cites[pmid])

# Clean dirty data
print("Cleaning metadata")
i = 0
for k in list(meta.keys()):
    i += 1
    if i % 1000000 == 0:
        print(i)
    try:
        meta[k]['doi'] = cleanString(meta[k]['doi'])
    except (KeyError):
        meta[k]['doi'] = strEmpty(None)
    try:
        temp = meta[k]['year']
    except (KeyError):
        meta[k]['year'] = strEmpty(None)
    try:
        meta[k]['journal'] = cleanString(meta[k]['journal'])
    except (KeyError):
        meta[k]['journal'] = strEmpty(None)
    try:
        meta[k]['issn'] = cleanString(meta[k]['issn'])
    except (KeyError):
        meta[k]['issn'] = strEmpty(None)
    try:
        temp = meta[k]['res']
    except (KeyError):
        meta[k]['res'] = strEmpty(None)
    try:
        temp = meta[k]['jcr']
    except (KeyError):
        meta[k]['jcr'] = strEmpty(None)
    try:
        temp = meta[k]['citation_count']
    except (KeyError):
        meta[k]['citation_count'] = 0
    try:
        temp = meta[k]['ref_count']
    except (KeyError):
        meta[k]['ref_count'] = 0

# Write metadata
# Metadata tsv
print("Writing metadata to tsv")
with open(metaOutPath + '.tsv', 'w') as file:
    for k in list(meta.keys()):
        outline = '\t'.join([str(k), meta[k]['doi'], str(meta[k]['year']), meta[k]['journal'], meta[k]['issn'], str(meta[k]['res']), str(meta[k]['jcr']), str(meta[k]['citation_count']), str(meta[k]['ref_count'])]) + '\n'
        file.write(outline)


# Citing network tsv
print("Writing citing data to tsv")
with open(citesOutPath + '.tsv', 'w') as file:
    for k in list(cites.keys()):
        outline = str(k) + '\t' + '|'.join([str(x) for x in cites[k]]) + '\n'
        file.write(outline)

# Cited network tsv
print("Writing citing data to tsv")
with open(citedByOutPath + '.tsv', 'w') as file:
    for k in list(citedBy.keys()):
        outline = str(k) + '\t' + '|'.join([str(x) for x in citedBy[k]]) + '\n'
        file.write(outline)


