start.year = 1995
# Get current indexed time from percent_indexed_in_icite.r
indexed.time = 2014 + read.table("./indexed_time.txt", header=F)[,1]
nih.papers = unique(read.csv("./publink 1980-2015 NIH only.csv", stringsAsFactors=F)[,"PMID"])
icite <- read.delim("./icite table.tsv", quote="", stringsAsFactors=FALSE, colClasses=c("numeric", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"), comment.char="", sep="\t")
nihflag = icite$year >= start.year & icite$pmid %in% nih.papers

quant=NULL
for (y in start.year:(floor(indexed.time))) {
  print(y)
  q = quantile(icite[icite$year == y & nihflag, "rcr"], c(0,0.1,0.2,0.4,0.5,0.6,0.8,0.9,0.95,0.99,0.999,1), na.rm=T)
  quant = rbind(quant, c(y, q))
}
q = quantile(icite[icite$year %in% start.year:(floor(indexed.time)) & nihflag, "rcr"], c(0,0.1,0.2,0.4,0.5,0.6,0.8,0.9,0.95,0.99,0.999,1), na.rm=T)
q2 = quantile(icite[icite$year %in% start.year:(floor(indexed.time-1)) & nihflag, "rcr"], c(0,0.1,0.2,0.4,0.5,0.6,0.8,0.9,0.95,0.99,0.999,1), na.rm=T)
quant = rbind(quant, q, q2)
q2
q
