meta = read.table("./mergeMeta.tsv", sep="\t", quote="", comment.char="",
                  colClasses=c("numeric", "character", "numeric", "character", "character",
                               "numeric", "numeric", "numeric", "numeric"), stringsAsFactors = F,
                  col.names=c('pmid', 'doi', 'year', 'journal', 'issn', 'res', 'jcr','citation_count', 'ref_count'))
meta75 = meta[meta$year >= 1975 & !(is.na(meta$year)),]
colnames(meta75)=c('pmid', 'doi', 'year', 'journal', 'issn', 'res', 'jcr','citation_count', 'ref_count')

# From start year until last complete year, read table of docs per year in dataset
start.year = 1995
end.year = 2014
x = data.frame(as.integer(names(tapply(meta75$year, meta75$year, length))), tapply(meta75$year, meta75$year, length))
y = x[x[,1] >= start.year & x[,1] <= end.year,]
# For sanity's sake, make x-axis start at 1
y[,1] = 1:nrow(y)

mylm = glm(y[,2] ~ y[,1], family=Gamma(link=log))

# Extrapolate to current year
z = exp(mylm$coefficients[1]) * exp(mylm$coefficients[2] * 1:(nrow(y) + 1))


percent.through.current.year = x[x[,1] == (end.year + 1), 2] / z[length(z)]
write.table(percent.through.current.year, "./indexed_time.txt", row.names=F, col.names=F)

