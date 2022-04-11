# the data are from here:
#https://www.scrna-tools.org/table

d = read.csv('input/scrna-tools.tableExport.csv')
dim(d)
table(d$Platform)
pls = strsplit(d$Platform,'/',fixed = T)
langs = sort(table(unlist(pls)),decreasing = T)
d$year = as.numeric(substr(d$Pub.Dates,1,4))
table(d$year,useNA = 'always')
table(substr(d$Added,1,4),useNA = 'always')
d$ayear = as.numeric(substr(d$Added,1,4))
table(d$year,d$ayear,useNA = 'always')


m = matrix(0,nrow=nrow(d),ncol=length(langs))
colnames(m) = names(langs)
for(i in 1:length(pls))
  m[i,pls[[i]]] = 1

mn = sweep(m,1,rowSums(m),'/')


stat = sapply(split.data.frame(mn,d$year),colSums)
stta = sapply(split.data.frame(mn,d$ayear),colSums)

statn = sweep(stat,2,colSums(stat),'/')
sttan = sweep(stta,2,colSums(stta),'/')

barplot(statn[1:10,],xlim=c(0,23),las=2,legend.text = T,col=rainbow(10),ylim=c(0,1))
barplot(sttan[1:10,],xlim=c(0,13),las=2,legend.text = T,col=rainbow(10),ylim=c(0,1))
