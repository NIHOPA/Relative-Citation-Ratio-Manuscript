library(quantreg)
library(jsonlite)
# Read current list of NIH-funded papers
start.year = 1980
# Input current indexed time here. = Last full year indexed plus fraction of current year
indexed.time = 2014 + read.table("./indexed_time.txt", header=F)[,1]
# Filter the last year if it's partial
acceptable.months = c("Jan", "Feb", "Mar", "Apr", "May", "Jun")
nih.papers = unique(read.csv("./publink 1980-2015 NIH only.csv", stringsAsFactors=F)[,"PMID"])
metadata = read.table("./mergeMeta.tsv", sep="\t", quote="", comment.char="",
                      colClasses=c("numeric", "character", "numeric", "character", "character",
                                   "numeric", "numeric", "numeric", "numeric"), stringsAsFactors = F, header=T,
                      col.names=c('pmid', 'doi', 'year', 'journal', 'issn', 'res', 'jcr','citation_count', 'ref_count'))

metadata = metadata[!is.na(metadata$pmid) & metadata$year >= start.year & !is.na(metadata$year),]

# Check year data
summary(metadata$year)

nihdata = metadata[metadata$year >= start.year & metadata$pmid %in% nih.papers,]
nihflag = metadata$year >= start.year & metadata$pmid %in% nih.papers

# Read FCRs in and match
# Use JCRs if FCRs aren't available (mostly uncited articles)
jcrmatch = read.csv("./jcrmap.csv", header = F, stringsAsFactors=F,
                    col.names = c("pmid", "jcr"))
fcrmatch = read.table("./fcrmap.tsv", sep="\t", header=F,
                      stringsAsFactors=F, col.names=c("pmid", "fcr"))
# Default FCR is NA
nihdata$fcr = rep(as.numeric(NA), nrow(nihdata))
summary(nihdata$fcr)
length(nihdata$fcr)
# Fill in JCR for FCR as secondary default
tempmatches = match(nihdata$pmid, jcrmatch$pmid)
nihdata$fcr[!is.na(tempmatches)] <- jcrmatch$jcr[na.omit(tempmatches)]
summary(nihdata$fcr)
# Fill in FCR for FCR if available
tempmatches = match(nihdata$pmid, fcrmatch$pmid)
nihdata$fcr[!is.na(tempmatches)] <- fcrmatch$fcr[na.omit(tempmatches)]
summary(nihdata$fcr)

# Check summary stats
tapply(nihdata$citation_count, nihdata$year, function(x) mean(x, na.rm=T))
tapply(nihdata$citation_count, nihdata$year, length)

# For month data, consult parsed Medline data
load("./medline15.meta.parsed.rdata")
lastyear = floor(indexed.time)
medline.sub = medline15.meta.parsed[medline15.meta.parsed$year == lastyear,]
# Extract month data
medline.sub$mo1 = substr(medline.sub$date, 5, 7)
medline.sub$mo2 = substr(medline.sub$date, 6, 8)
acceptable.lastyear.pmids = medline.sub[medline.sub$mo1 %in% acceptable.months | medline.sub$mo2 %in% acceptable.months, "pmid"]
# acceptable.lastyear.pmids = metadata[metadata$year == lastyear, "pmid"]
# Generate regressions
regs = NULL
year.pmids = list()
for (y in start.year:(floor(indexed.time))) {
  # Subset & remove entries with missing data
  temp = nihdata[nihdata$year == y,]
  temp = temp[!is.na(temp$pmid) & !is.na(temp$citation_count) & !is.na(temp$fcr),]
  # If we've restricted the month range for the last year, filter here
  if (y == lastyear) temp = temp[temp$pmid %in% acceptable.lastyear.pmids,]
  cpy = temp$citation_count / (indexed.time - y)
  print(paste("Year", y, "documents", nrow(temp)))
  # Don't time-normalize for partial year at the end; you'll get whacky values
  if (y == lastyear) cpy = temp$citation_count
  fcr = temp$fcr
  # Quantile regression
  this.rq = rq(cpy ~ fcr, R=diag(2), r=c(0,0), method="fnc")
  regs = rbind(regs, data.frame("Year"=y, "Intercept"=this.rq$coef[1], "Slope"=this.rq$coef[2]))
  # Median should be very close to 1.0
  print(median(cpy / (this.rq$coef[1] + this.rq$coef[2] * fcr)))
  year.pmids[[as.character(y)]] = temp$pmid
}
write(toJSON(year.pmids), file="./pmids_used.json")
write.csv(regs, "./qregs update.csv", row.names=F)

# Fill in ECRs for NIH data
nihdata$ecr = rep(as.numeric(NA), nrow(nihdata))
tempmatches = match(nihdata$year, regs$Year)
nihdata$ecr[!is.na(tempmatches)] <- regs$Intercept[na.omit(tempmatches)] + (regs$Slope[na.omit(tempmatches)] * nihdata$fcr[!is.na(tempmatches)])
summary(nihdata$ecr)

# Calculate NIH CPYs
cpy.denom = indexed.time - nihdata$year
summary(cpy.denom)
cpy.denom[which(cpy.denom < 1)] = 1
summary(cpy.denom)
nihdata$cpy = nihdata$citation_count / cpy.denom
summary(nihdata$cpy)

# Calculate NIH RCRs
nihdata$rcr = nihdata$cpy / nihdata$ecr
summary(nihdata$rcr)

# Fill in metadata with FCR, ECR, RCR
# Default FCR is NA
metadata$fcr = rep(as.numeric(NA), nrow(metadata))
summary(metadata$fcr)
length(metadata$fcr)
# Fill in JCR for FCR as secondary default
tempmatches = match(metadata$pmid, jcrmatch$pmid)
metadata$fcr[!is.na(tempmatches)] <- jcrmatch$jcr[na.omit(tempmatches)]
summary(metadata$fcr)
# Fill in FCR for FCR if available
tempmatches = match(metadata$pmid, fcrmatch$pmid)
metadata$fcr[!is.na(tempmatches)] <- fcrmatch$fcr[na.omit(tempmatches)]
summary(metadata$fcr)

# Calculate ECRs
metadata$ecr = rep(as.numeric(NA), nrow(metadata))
tempmatches = match(metadata$year, regs$Year)
metadata$ecr[!is.na(tempmatches)] <- regs$Intercept[na.omit(tempmatches)] + (regs$Slope[na.omit(tempmatches)] * metadata$fcr[!is.na(tempmatches)])
summary(metadata$ecr)

# Calculate CPYs
cpy.denom = indexed.time - metadata$year
summary(cpy.denom)
cpy.denom[which(cpy.denom < 1)] = 1
summary(cpy.denom)
metadata$cpy = metadata$citation_count / cpy.denom
summary(metadata$cpy)

# Calculate RCRs
metadata$rcr = metadata$cpy / metadata$ecr
summary(metadata$rcr)
quantile(metadata$rcr, 1:199/200, na.rm=T)

# Truncate final year if partial year data
flag = metadata$year == lastyear & !(metadata$pmid %in% acceptable.lastyear.pmids)
metadata$fcr[flag] = NA
metadata$ecr[flag] = NA
metadata$rcr[flag] = NA
flag = nihdata$year == lastyear & !(nihdata$pmid %in% acceptable.lastyear.pmids)
nihdata$fcr[flag] = NA
nihdata$ecr[flag] = NA
nihdata$rcr[flag] = NA


# Uncited papers have RCR = 0
flag = metadata$citation_count == 0 & metadata$year <= lastyear & (metadata$year < lastyear | (metadata$year == lastyear & metadata$pmid %in% acceptable.lastyear.pmids))
metadata$rcr[flag] = 0
flag = nihdata$citation_count == 0 & nihdata$year <= lastyear & (nihdata$year < lastyear | (nihdata$year == lastyear & nihdata$pmid %in% acceptable.lastyear.pmids))
nihdata$rcr[flag] = 0
summary(metadata$rcr)
summary(metadata$ecr)
summary(metadata$fcr)
summary(metadata$cpy)
summary(nihdata$rcr)
summary(nihdata$ecr)
summary(nihdata$fcr)
summary(nihdata$cpy)

summary(metadata$rcr[nihflag])
summary(metadata$ecr[nihflag])
summary(metadata$fcr[nihflag])
summary(metadata$cpy[nihflag])

# Save
write.table(metadata, "./icite table.tsv", sep="\t", quote=F,
            na="", row.names=F, col.names=T)
