library(stringr)
library(gh)

# Limitations
# 1. I didn't check whether github link correspond to the toold developed for this publication, probably it is just tool used in work

# check md5sums #####
e = read.table('processed/tars/sum')
o = read.table('processed/tars/check.sums')
i = which(o$V1 != e$V1[match(o$V2,e$V2)])
o[i,]

# parse zgrep results ####
sids = sub('.tar.gz','',list.files('processed/tars',pattern = '*gz'))
done = sub('.github.txt','',list.files('processed/ghmatches/'))
done
setdiff(sids,gsub('.gz','',done))

# files are in json, but extension is xml...
# depracated: files are too large to load them to memory
# pmc2gh = lapply(done,function(f){
#   t = readLines(paste0('processed/ghmatches/',f,'.github.txt'),skipNul = T)
#   pmcids = stringr::str_match(t,'PMC\\d+.xml')
#   ghlinks = stringr::str_match_all(t,'[^\\s()"\']*github.com[^\\s()"\']*')
#   ghlinks = lapply(ghlinks,function(x)sub('\\.$|,$','',x))
#   names(ghlinks) = pmcids
#   ghlinks
# })
# names(pmc2gh) = sub('_json_unicode','',done)
# sapply(pmc2gh,length)

# prepare for individual file extraction 
# deprecated: no need for it, I'll use zgrep output
# dir.create('processed/jsons')
# for(n in names(pmc2gh)){
#   path = paste0('processed/jsons/',n)
#   if(!dir.exists(path)) dir.create(path)
#   writeLines(names(pmc2gh[[n]]),paste0(path,'/files.txt'))
# }


# 
parseTarExtracts = function(f){
  input = file(f, "r")
  r = list()
  repeat{
    l = readLines(input, n = 1,skipNul = TRUE)
    if(is.null(l) || length(l) == 0 || is.na(l) || l =='' ) break;
    pmc = substr(l,1,14)
    ghlinks = str_match_all(l,'[^\\s()"\']*github.com[^\\s()"\']*')[[1]][,1]
    ghlinks = sub('\\.$|,$','',ghlinks)
    l = gsub('^[^{]+|[^}]*$','',l)
    js = jsonlite::parse_json(l)
    art = js$documents[[1]]$passages[[1]]
    art$ghlinks = ghlinks
    r[[pmc]] = art
    cat('\r',length(r),'     ')
  }
  cat('\n')
  close(input)
  r
}

fls = list.files('processed/ghmatches/')
arts = lapply(paste0('processed/ghmatches/',fls),parseTarExtracts)
sum(sapply(arts,length))


#saveRDS(arts,'rds/arts.rds')
arts = readRDS('rds/arts.rds')
arts = do.call(c,arts)
length(arts)

# keyywords
kwds = sapply(arts,function(a){r=a$infons$kwd;if(is.null(r)){r=NA};r})
table(is.na(kwds))
t = sort(table(unlist(strsplit(kwds,' '))),dec=TRUE)
t[1:10]
par(mar=c(3,10,1,1),cex=0.5)
barplot(t[40:1],horiz = T,las=2)
length(t)
table(t>500)

# I used just "github" for zgrep (without .com)
l = sapply(arts,function(x)length(x$ghlinks))
table(l)
arts = arts[l>0]

years = sapply(arts,function(x)x$infons$year)
table(sapply(years,length))
# 0     1 
# 3 75622 
arts = arts[sapply(years,length)==1] 

years = as.numeric(sapply(arts,function(x)x$infons$year))
range(years)
hist(years,2009:2022)

names(arts) = sub('.xml','',names(arts))

ghs = lapply(arts, function(x){
  t = sapply(strsplit(x$ghlinks,'github.com/+'),'[',2)
  r = unique(as.data.frame(do.call(rbind,lapply(strsplit(t,'/+'),'[',1:2))))
  colnames(r) = c('user','rep')
  r = r[!is.na(r$rep) & !is.na(r$user),]
  r$rep = gsub('\\.git|]|\\\\n|;','',r$rep)
  r$rep = gsub('\\\\|<|\\\\t','',r$rep)
  unique(r)
  })

table(sapply(ghs,nrow))
f = sapply(ghs,nrow) > 0
arts = arts[f]
ghs = ghs[f]
for(n in names(ghs)) ghs[[n]]$pmc = n

ghs = do.call(rbind,ghs)
ghs$id = paste0(ghs$user,'.',ghs$rep)
years = as.numeric(sapply(arts,function(x)x$infons$year))
names(years) = names(arts)
ghs$year = years[ghs$pmc]

ghs[1:2,]
table(table(ghs$id))
#saveRDS(ghs,'rds/gh2pmc.rds')

ghs = readRDS('rds/gh2pmc.rds')
# chose earliest repo mention
ughs = unique(ghs[,c('id','user','rep')])
dim(ughs)

# use github token here
token=''
langs = list()
if(file.exists('rds/langs.temp.rds'))
  langs = readRDS('rds/langs.temp.rds')

qn = 0
for(i in 1:nrow(ughs)){
  if(length(langs) >= i && !is.null(langs[[i]]) && !is.null(langs[[i]]$langs) && !is.null(langs[[i]]$info)) next
  # not more than 5000 request per hour are allowed
  # I'll wait before start to make sure, that results are not affected by previous runs (remove "i!=1")
  # I'll wait one hour, to make sure that i'm not blocked
  #if(i!= 1 && i %% 2450 == 1){
  if(qn == 2450){
    t = Sys.time()
    repeat{
      wt = 3610 - ceiling(as.double(difftime(Sys.time(),t,units = 'secs')))
      cat('\r',i,' wait for ',wt,' secs     ')
      Sys.sleep(10)
      if(wt < 0) break
    }
    qn = 0
  }
  cat('\r',i,qn,'    ')
  r = d = NULL
  qn = qn + 1
  try({r=gh(paste0("GET http://api.github.com/repos/",ughs$user[i],"/",ughs$rep[i],"/languages"),.token=token)},silent = TRUE)
  try({d=gh(paste0("GET http://api.github.com/repos/",ughs$user[i],"/",ughs$rep[i]),.token=token)},silent = TRUE)
  langs[[i]] = list(langs=r,info=d)
  if(i %% 1000 == 0)
    saveRDS(langs,'rds/langs.temp.rds')
}
names(langs) = ughs$id
saveRDS(langs,'rds/langs.filtered.rds')

# saveRDS(arts,'rds/arts.filtered.rds')

# filter #####
ghlist = readRDS('rds/langs.rds')
#langs = readRDS('rds/langs.temp.rds')
gh2art = readRDS('rds/gh2pmc.rds')
artslist = readRDS('rds/arts.filtered.rds')

# filter to retain only real repos
# add number of mentioned repos to articles: to make it possible to filter by this parameter
length(artslist)
length(ghlist)
nrow(gh2art)

# check repos that do not have data
inx=which(sapply(ghlist,function(x)is.null(x$langs)))
plot(inx,t='l')
# more or less linear, so likely there were no problems with rate limit
length(inx)


# make lang usage matrix
ghlist = ghlist[-inx]
lnames = sort(table(unlist(lapply(ghlist,function(x)names(x$langs)))),dec=T)
length(lnames)
lnames[1:10]


lang.mat = matrix(0,ncol=length(lnames),nrow=length(ghlist),dimnames=list(names(ghlist),names(lnames)))
for(i in 1:length(ghlist)){
  try({lang.mat[i,names(ghlist[[i]]$langs)] = as.numeric(ghlist[[i]]$langs)},silent=TRUE)
}
lang.mat[1:5,1:5]

# repos with no langs
inxs = which(apply(lang.mat,1,sum)==0)
length(inxs)
# probably need to check how is it possible


lang.mat = lang.mat[-inxs,]
ghlist = ghlist[-inxs]
dim(lang.mat)
length(ghlist)

flds = c('created_at','language','updated_at','pushed_at','git_url')
gh = lapply(ghlist,function(g){as.data.frame(lapply(g$info[flds],function(x)ifelse(is.null(x),NA,x)))})
gh = do.call(rbind,gh)
dim(gh)
apply(is.na(gh),2,sum)
sort(table(gh$language))

gh$most.freq.lang = colnames(lang.mat)[apply(lang.mat,1,which.max)]
table(gh$most.freq.lang == gh$language)
gh[1:2,]
gh$art.count =  as.integer(table(gh2art$id)[rownames(gh)])
table(gh$art.count,useNA = 'always')
gh[order(gh$art.count,decreasing = T)[1:10],]

table(gh2art$id %in% rownames(gh))
gh2art = gh2art[gh2art$id %in% rownames(gh),]

# find earliest paper
eart = do.call(rbind,lapply(split(gh2art,gh2art$id),function(x)x[order(x$year)[1],]))
eart = eart[rownames(gh),]
gh$earliest.pmc = eart$pmc
gh$earliest.pmc.year = eart$year
table(gh$earliest.pmc.year,substr(gh$created_at,1,4))

# prepare article table
f = function(x)ifelse(is.null(x),NA,x)
arts = do.call(rbind,lapply(artslist,function(a)
  data.frame(pmid=f(a$infons$`article-id_pmid`),
             year=f(a$infons$year),
             title=f(a$text),
             keywords=f(a$infons$kwd))
))
arts[1:2,]
dim(arts)
table(rownames(arts) %in% gh2art$pmc)
arts = arts[rownames(arts) %in% gh2art$pmc,]
arts[1:2,]
arts$gh.count = as.numeric(table(gh2art$pmc)[rownames(arts)])
table(arts$gh.count)

# save
# saveRDS(arts,'rds/arts.table.final.rds')
# saveRDS(gh2art,'rds/gh2art.final.rds')
# saveRDS(gh,'rds/gh.table.final.rds')
# saveRDS(lang.mat,'rds/lang.mat.rds')
