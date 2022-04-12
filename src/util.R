#par(mfrow=c(2,2),tcl=-0.2,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1),bty='n')
#source("/home/mazin/R/util.R")
#yaxs='i'
#library(RMySQL)
#con = dbConnect(MySQL(), user="reader",dbname="species75jl", host="127.0.0.1")

#debug();findLineNum();setBreakpoint();browser()

`%c%` = function(x,y)paste0(x,y)

log10p1 = function(x)log1p(x)/log(10)

calcPSI = function(i,e,l,lr,min.cov=10){
  f = i + e < min.cov
  i = i/(l+lr-1)
  e = e/(lr-1)
  psi = i/(i+e)
  psi[f] = NA
  psi
}



readNamedMM = function(f){
  require(Matrix)
  if(file.exists(paste0(f,'.mtx'))){
    m = Matrix::readMM(paste0(f,'.mtx'))
    rownames(m) = readLines(paste0(f,'_key1.csv'))
    colnames(m) = readLines(paste0(f,'_key2.csv'))
  }else{
    m = Matrix::readMM(paste0(f,'.mtx.gz'))
    rownames(m) = readLines(paste0(f,'_key1.csv.gz'))
    colnames(m) = readLines(paste0(f,'_key2.csv.gz'))
  }
  m
}

pcor=function(m1,m2) crossprod(scale(m1), scale(m2))/(nrow(m1)-1)

plot.as.hm = function(x,y,xbins=100,ybins=100,cols=c('white','gray','blue','orange','red'),zfun=identity,leg.title='',num.leg.tic=NULL,legend=TRUE,trimZq=0,xlim=NULL,ylim=NULL,new=TRUE,xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),...){
  xlab;ylab;
  f = !is.na(x) & !is.infinite(x) & !is.na(y) & !is.infinite(y) 
  if(!is.null(xlim)) f = f & x >= xlim[1] & x <= xlim[2]
  if(!is.null(ylim)) f = f & y >= ylim[1] & y <= ylim[2]
  x = x[f]
  y = y[f]
  if(length(xbins)==1) xbins = seq(min(x),max(x),length.out = xbins+1)
  if(length(ybins)==1) ybins = seq(min(y),max(y),length.out = ybins+1)
  
  f = x >= xbins[1] & x <= xbins[length(xbins)] & y >= ybins[1] & y <= ybins[length(ybins)]
  x = x[f]
  y = y[f]
  
  xb = findInterval(x,xbins,rightmost.closed = TRUE)
  yb = findInterval(y,ybins,rightmost.closed = TRUE)
  
  
  z = as.matrix(table(factor(yb,levels=1:(length(ybins)-1)),factor(xb,levels=1:(length(xbins)-1))))
  z = trimQ(z,trimZq)
  if(is.null(xlim)) xlim = range(x)
  if(is.null(ylim)) ylim = range(y)
  if(new)
    plot(1,t='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
  z2col=function(x)num2col(zfun(x),cols)
  rect(rep(xbins[-length(xbins)],each =length(ybins)-1),
       rep(ybins[-length(ybins)],times=length(xbins)-1),
       rep(xbins[-1]            ,each =length(ybins)-1),
       rep(ybins[-1]            ,times=length(xbins)-1),
       col=z2col(z),border = NA)
  if(legend)
    plotColorLegend2(grconvertX(1,'npc','nfc'),1,grconvertY(0,'npc','nfc'),grconvertY(1,'npc','nfc'),
                   range(z),range(z),zfun,z2col=z2col,leg=num.leg.tic,title=leg.title)
  invisible(list(xbins=xbins,ybins=ybins,z=z))
}

trimQ = function(x,q){
  if(q==0) return(x)
  qq=quantile(x,sort(c(q,1-q)))
  x[x<=qq[1]] = qq[1]
  x[x>=qq[2]] = qq[2]
  x
}


castXYtable = function(x,y,i){
  ys = sort(unique(y))
  xs = sort(unique(x))
  m = matrix(NA,nrow=length(xs),ncol=length(ys),dimnames=list(xs,ys))
  x = as.character(x)
  y = as.character(y)
  for(j in 1:length(x))
    m[x[j],y[j]]= i[j]
  m
}

plotVolcano = function(l2fc,pv=NULL,gnames=NULL,pch=19,cex=c(0.5,2),xlab='l2FC',bty='n',fdr.thr=0.05,topN2name=10,topHigherThan='pv',glabcol='red',plot.pv.thr=TRUE,gncex=1,gnfont=2,...){
  if(is.null(pv)){
    pv = l2fc$pv
    gnames = rownames(l2fc)
    l2fc = l2fc$l2fc
  }
  fdr = p.adjust(pv,m='BH')
  if(!is.null(topHigherThan) & topHigherThan=='pv')
    topN2name = min(topN2name,sum(pv<0.05))
  if(!is.null(topHigherThan) & topHigherThan=='fdr')
    topN2name = min(topN2name,sum(fdr<0.05))
  pv = log10(pv)
  fdr = max(pv[fdr<fdr.thr])
  
  pvcex = -pv
  pvcex = (pvcex-min(pvcex))/(max(pvcex)-min(pvcex))*(cex[2]-cex[1]) + cex[1]
  plot(l2fc,pv,bty=bty,pch=pch,ylim=rev(range(pv,fdr,log10(0.05))),xlab=xlab,ylab='log10(pv)',cex=pvcex,...)
  
  inx = order(pv)[1:topN2name]
  par(xpd=NA)
  text(l2fc[inx],pv[inx],gnames[inx],col=glabcol,adj=c(-0.1,-0.1),cex=gncex[inx],font=gnfont[inx])
  par(xpd=FALSE)
  if(plot.pv.thr){
    hor = c(log10(0.05),fdr)
    txt = c('P = 0.05',paste0('FDR = ',fdr.thr))
  }else{
    hor = c(fdr)
    txt = paste0('FDR = ',fdr.thr)
  }
  abline(h=hor,col='gray',lty=2)
  text(grconvertX(1,'npc','user'),hor,txt,adj=c(1,1.1))
  abline(v=0,col='gray',lty=2)
}

char2col = function(t,bpal='Set1',colfun=rainbow,palette=TRUE){
  torig = t
  t = sort(unique(t))
  suppressWarnings({
    if(all(!is.na(as.integer(t)))){
      t = as.character(sort(as.integer(t)))
    }
  })
  if(length(t)<10)
    r=setNames(RColorBrewer::brewer.pal(max(3,length(t)),bpal),t)[1:length(t)]
  else
    r=setNames(colfun(length(t)),t)
  if(!palette)
    r = r[torig]
  r
}

recycle = function(v,i){
	v[(i-1)%%length(v)+1]
}

addExonNumber = function(segs){
	print('\taddExonNumber')
	require(data.table)
	ss = split(segs,segs$gene_id)
	k = 1
	r = lapply(ss,function(g){
		cat('\r',k,length(ss))
		k <<- k + 1
		g = g[order(g$start*g$strand),]
		g$exon.number = NA
		n = 1
		exn = TRUE
		for(i in 1:nrow(g)){
			if(exn){
				if(substr(g$sites[i],2,2) %in% c('d','.'))
					exn = FALSE
			}else{
				if(substr(g$sites[i],1,1) %in% c('a','.')){
					exn = TRUE
					n = n + 1
					if(substr(g$sites[i],2,2) %in% c('d','.'))
						exn = FALSE
				}
			}
			g$exon.number[i] = n
			if(substr(g$sites[i],1,1) %in% c('a','.') & g$position[i] == 'LAST'){ # to account for internal lasts
				n = n - 1
				exn = FALSE
			}
		}
		g$seg.id = rownames(g)
		g
	})
	r = data.frame(rbindlist(r))
	rownames(r) = r$seg.id
	r$seg.id = NULL
	r[rownames(segs),]
}

num2col = function(d,col=c('blue','cyan','gray','orange','red'),minx=min(d),maxx=max(d)){
  if(sd(d)==0)
    return(rep(col[1],length(d)))
	d[d<minx] = minx
	d[d>maxx] = maxx
	bwr = colorRamp(col,alpha=TRUE)
	apply(bwr(scaleTo(d,0,1,minx=minx,maxx=maxx))/255,1,function(x)rgb(x[1],x[2],x[3],x[4]))
}

isColors <- function(x) {
  x %in% colors() | (substr(x,1,1)=='#' & (nchar(x)==7 |  nchar(x)==9))
  # sapply(x, function(X) {
  #   tryCatch(is.matrix(col2rgb(X)), 
  #            error = function(e) FALSE)
  # })
}

imageTable = function(x,y=NULL,z,z2col=num2col,zfun=identity,zlim=NULL,step=1,ystep=step,brdcol=NA,legN=5,ylim=NULL,xlim=NULL,plot.legend=FALSE,legend.args=list(),num.leg.tic=NULL,...){
  if(is.null(y)){
    y = x$y
    x = x$x
  }
  # make color
  if(is.numeric(z)){
    if(!is.null(zlim)){
      z[z>zlim[2]] = zlim[2]
      z[z<zlim[1]] = zlim[1]
    }else
      zlim=range(z)
    zorig = z
    z = zfun(z)
    col = z2col(c(zfun(zlim),z))[-(1:length(zlim))]
  }else if(all(isColors(z))){
    plot.legend = FALSE
    col = z
  }else{
    if(!is.character(z2col)){
      z2col = char2col(z)
    }
    col = z2col[z]
  }
  if(is.null(xlim)) xlim = c(min(x)-step/2,max(x)+step/2)
  if(is.null(ylim)) ylim = c(min(y)-ystep/2,max(y)+ystep/2)
  plot(1,t='n',xlim=xlim,ylim=ylim,...)
  rect(x-step/2,y-ystep/2,x+step/2,y+ystep/2,col=col,border=brdcol)
  if(plot.legend){
    if(is.character(z2col)){ 
      # categorical
      if(is.null(legend.args$x)){
        legend.args$x = grconvertX(1,'npc','user')
        legend.args$y = grconvertY(1,'npc','user')
      }
      legend.args.def = list(xpd=NA,pch=19,col=z2col,legend=names(z2col),bty=par('bty'))
      for(n in names(legend.args.def))
        if(is.null(legend.args[[n]])) legend.args[[n]] = legend.args.def[[n]]
      do.call(legend,legend.args)
    }else{
      # numerical
      plotColorLegend2(grconvertX(1,'npc','nfc'),1,grconvertY(0,'npc','nfc'),grconvertY(1,'npc','nfc'),zlim,range(zorig),zfun,z2col,leg=num.leg.tic,title=legend.args$title)
    }
  }
}

imageTable_old = function(x,y=NULL,col,step=1,brdcol=NA,colgrad=c('blue','cyan','gray','orange','red'),col.fun=identity,plot.legend=FALSE,legN=5,ylim=NULL,xlim=NULL,...){
	if(is.null(y)){
		y = x$y
		x = x$x
	}
	col.orig = col
	col = col.fun(col)
	if(!is.character(col) && !is.integer(col) && any(col!=floor(col))){
		col = num2col(col,colgrad)
	}
	if(is.null(xlim)) xlim = c(min(x)-step/2,max(x)+step/2)
	if(is.null(ylim)) ylim = c(min(y)-step/2,max(y)+step/2)
	plot(1,t='n',xlim=xlim,ylim=ylim,...)
	rect(x-step/2,y-step/2,x+step/2,y+step/2,col=col,border=brdcol)
	if(plot.legend){
		legN = legN - 1
		cr = range(col.orig,na.rm=T)
		crt = col.fun(cr)
		
		lcols = seq(crt[1],crt[2],length.out=200)
		at = seq(cr[1],cr[2],length.out=20000)
		att=ceiling(scaleTo(col.fun(at))*legN)
		leg = at[c(1,sapply(1:legN,function(i)max(which(att==i))))]
		leg = leg[!is.na(leg)]
		leg[leg>10] = round(leg[leg>10])
		at = findInterval(col.fun(leg),lcols,all.inside = T)
		if(max(leg)/min(leg)>100)
			leg = gsub('(e[+-])0+','\\1',format(leg,sciencific=T,digits=1,trim=T))
		leg[leg=='0e+' | leg=='0e-'] = '0'
		plotColorLegend(grconvertX(max(x)+step/2,'user','nfc'),1,0.4,0.8, num2col(lcols,colgrad),at,leg)
	}
}

imageWithText = function(d,t=NULL,digits=2,text.col=NULL,xaxlab=rownames(d),yaxlab=colnames(d),centerColors0=TRUE,...){
	if(is.null(t))
		t = round(d,digits = digits)
	pars = list(...)
	if(is.null(pars$col)) pars$col= getPal()
	pars$x = 1:nrow(d)
	pars$y = 1:ncol(d)
	pars$z = d
	pars$xaxt='n'
	pars$yaxt='n'
	if(!is.null(col) & centerColors0){
	  m = max(abs(d),na.rm=T)
	  if(is.numeric(centerColors0))
	    m = centerColors0
	  pars$breaks = seq(-m,m,length.out = length(pars$col)+1)
	}
	do.call(image,pars)
	if(is.null(text.col))
		text.col = 'black'
	text(rep(pars$x,times=length(pars$y)),rep(pars$y,each=length(pars$x)),t,col=text.col)
	if(!is.null(xaxlab))
		axis(1,pars$x,xaxlab)
	if(!is.null(yaxlab))
		axis(2,pars$y,yaxlab)
}

plotTranscripts = function(a,
													 ylim=c(0,length(unique(a$transcript_id))),
													 xlim=c(ifelse(a$strand[1]=='+',min(a$start),max(a$stop)),ifelse(a$strand[1]=='+',max(a$stop),min(a$start))),
													 xlab=a$chr_id[1],
													 new=TRUE,yspace=0.8,exon.col='black',cds.col='black',...){
	transc = split(a,a$transcript_id)
	transc = transc[order(sapply(transc,function(x){max(x$stop)-min(x$start)}))]
	if(new)
		plot(1,t='n',xlim=xlim,ylim=ylim,yaxt='n',ylab='',xlab=xlab,...)
	ystep = (ylim[2]-ylim[1])/length(transc)
	for(i in 1:length(transc)){
		y = ylim[1] + ystep*i - ystep/2
		t = transc[[i]]
		lines(c(min(t$start),max(t$stop)),c(y,y))
		f = t$feature == 'exon'
		rect(t$start[f],y-ystep/2*yspace,t$stop[f],y+ystep/2*yspace,col = 'white',border = exon.col)
		f = t$feature == 'CDS'
		if(sum(f)>0)
			rect(t$start[f],y-ystep/2*yspace,t$stop[f],y+ystep/2*yspace,col = cds.col,border=cds.col)
	}
	text(par('usr')[1],seq(ylim[1]+ystep/2,by = ystep,length.out = length(transc)),sapply(transc,function(x)x$transcript_name[1]),adj = c(1,0.5),xpd=T)
}


loadEnsGTF = function(f,features=NULL){
	r = read.table(f,sep='\t')
	if(!is.null(features ))
		r = r[r$V3 %in% features,]
	a = lapply(strsplit(r$V9,';\\s?',perl=T),function(x){x=strsplit(x,'[ =]');setNames(sapply(x,'[',2),sapply(x,'[',1))})
	names = unique(unlist(lapply(a,names)))
	a = do.call(rbind,lapply(a,'[',names))
	colnames(a) = names
	r = r[,c(1:5,7)]
	colnames(r) = c('chr_id','type','feature','start','stop','strand')
	cbind(r,a)
}


parseHisat2LogP = function(f){
	if(!file.exists(f)) return(NULL)
	.f = function(l,p)
		as.numeric(gsub('%','',strsplit(l[grep(p,l,fixed = T)],' ')[[1]][1]))
	l = readLines(f)
	l = gsub('^ +','',l)
	fields = c(total='reads; of these:',
						 concord.once='aligned concordantly exactly 1 time',
						 concord.multi='aligned concordantly >1 times',
						 discordantly.once='aligned discordantly 1 time',
						 single.mate.once='aligned exactly 1 time',
						 single.mate.multi='aligned >1 times',
						 alignment.rate='overall alignment rate')
	sapply(fields,function(f).f(l,f))
}


parseHisat2LogS = function(f){
	if(!file.exists(f)) return(NULL)
	.f = function(l,p)
		as.numeric(gsub('%','',strsplit(l[grep(p,l,fixed = T)],' ')[[1]][1]))
	l = readLines(f)
	l = gsub('^ +','',l)
	fields = c(total='reads; of these:',
						 once='aligned exactly 1 time',
						 multi='aligned >1 times',
						 alignment.rate='overall alignment rate')
	sapply(fields,function(f).f(l,f))
}

makeTransparent = function(..., alpha=0.5) {
	
	if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
	
	alpha = floor(255*alpha)  
	newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
	
	.makeTransparent = function(col, alpha) {
		rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
	}
	
	newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
	
	return(newColor)
	
}

normValue = function(x,from=c(min(x,na.rm=T),max(x,na.rm=T)),to=c(0,1)){
	(x-from[1])/(from[2]-from[1])*(to[2]-to[1]) + to[1]
}

matches = function(x,y){
	x = data.frame(xid=1:length(x),v=x)
	y = data.frame(yid=1:length(y),v=y)
	r = merge(x,y,all.x=T)
	split(r$yid,r$xid)
}

suffix.matches = function(x,y,min.len=5){
	r = vector('list',length(x))
	ync = nchar(y)
	for(l in min.len:max(nchar(x))){
		f = nchar(x) == l
		r[f] = matches(x[f],substr(y,ync-l+1,ync))
	}
	r
}


plotPNG = function(f,x,y,w){
	# coordinates (x,y,w) are relative, x and y mark center of image, w is width, figures is not stratched!
  if(endsWith(tolower(f),'png')){
	  require(png)
	  img = readPNG(f)
  }
  if(endsWith(tolower(f),'jpg') | endsWith(tolower(f),'jpeg')){
    require(jpeg)
    img = readJPEG(f)
  }
	fy = grconvertY(0:1,'inches','npc')
	fx = grconvertX(0:1,'inches','npc')
	f = (fy[1]-fy[2])/(fx[1]-fx[2])
	h = w/dim(img)[2]*dim(img)[1]*f
	rasterImage(img,grconvertX(x-w/2,'npc','user'),grconvertY(y-h/2,'npc','user'),
							grconvertX(x+w/2,'npc','user'),grconvertY(y+h/2,'npc','user'))
}

my.binom.test = function(s,f=NULL){
	if(length(s)>1){
		f = s[2]
		s = s[1]
	}
	if(f+s == 0)
		return(c(NA,NA,NA))
	r=binom.test(s,s+f)
	c(r$estimate,r$conf.int)
}

revList = function(l){
	s = sapply(l,length)
	n = lapply(1:length(l),function(i){rep(names(l)[i],s[i])})
	split(unlist(n),unlist(l))
}

plot2Sign = function(s1,s2,v1,v2,n1,n2,main='',xlab=n1,ylab=n2,cols=c('gray',"#377EB8","#4DAF4A",'#E41A1C'),cex=rep(1,length(s1)),pch=19,range=NULL,ablineb=1,...){
	s1[is.na(s1)] = FALSE
	s2[is.na(s2)] = FALSE
	fa = function(x,y,f,cex,...){
		x = x[f]
		y = y[f]
		points(x,y,cex=cex[f],...)
		if(length(x)>1)
		  paste0("; r=",round(cor(x,y,use='pair'),2),'; N=',length(x),'; ',round(sum(sign(x) == sign(y),na.rm=T)/sum(!is.na(x-y))*100,0),'%',sep='')
		else
		  paste0("; r=NA; N=",length(x),'; ',round(sum(sign(x) == sign(y),na.rm=T)/sum(!is.na(x-y))*100,0),'%',sep='')
	}
	if(is.null(range))
		range = max(abs(range(v1,v2,na.rm=T)))
	range = c(-range,range)
	plot(1,t='n',xlim=range,ylim=range,xlab=xlab,ylab=ylab,main=main,...)
	cr = c(fa(v1,v2,!s1 & !s2,cex=cex,pch=pch,col=cols[1]),
				 fa(v1,v2, s1 & !s2,cex=cex,pch=pch,col=cols[2]),
				 fa(v1,v2,!s1 &  s2,cex=cex,pch=pch,col=cols[3]),
				 fa(v1,v2, s1 &  s2,cex=cex,pch=pch,col=cols[4]))
	abcol = '#666666'
	abline(v=0,lty=2,col=abcol)
	abline(h=0,lty=2,col=abcol)
	abline(a=0,b=ablineb,lty=2,col=abcol)
	#leg.names  = paste(c('!','!','',''),n1,' ',c('!','','!',''),n2,sep='')
	#leg.names = c('not sign.',paste('',n2),paste('',n1),'both')
	leg.names=c('none',n1,n2,'both')
	legend('topleft',col=cols[4:1],pch=pch,legend=paste(leg.names,cr,sep='')[4:1],bty='n')
	#legend('bottomright',col=cols[3:4],pch=pch,legend=paste(leg.names,cr,sep='')[3:4],bty='n')
	
}

first2Upper = function(t){
	paste0(toupper(substr(t,1,1)),substr(t,2,nchar(t)))
}

my.blast = function(db,seq,blast.arg = NULL,tmp.fa='~/tmp/blast/tmp.fa'){
	require(rBLAST)
	if(!is.null(db)){
		write.fasta(db,file=tmp.fa)
		makeblastdb(tmp.fa)
	}
	bl = blast(db=tmp.fa)
	if(!is.null(blast.arg)) blast.arg = paste0('-',names(blast.arg),' ',blast.arg,collapse = ' ')
	predict(bl, DNAStringSet(seq),BLAST_args = blast.arg)
}

plotColorLegend2 = function(x0,x1,y0,y1,fullzlim,zlim,zfun,z2col,N=100,ntic=5,leg=NULL,title=NULL){
  # make tics
  if(is.null(leg)){
    ztic = seq(zlim[1],zlim[2],length.out = 1e5)
    ztict = zfun(ztic)
    ztticat = seq(zfun(zlim[1]),zfun(zlim[2]),length.out = ntic)
    leg = ztic[findInterval(ztticat,ztict)]
    leg[1] = zlim[1]
    leg[ntic] = zlim[2]
  
    # adjust tic step
    d = (10^floor(log10(leg[2]-leg[1])))/2
    leg = unique(d*round(leg / d))
    leg = leg[leg>=zlim[1]]
  }
  # make colors
  ztat = zfun(leg)
  ztcol = sort(unique(c(ztat,seq(zfun(zlim[1]),zfun(zlim[2]),length.out = N))))
  col = z2col(c(zfun(fullzlim),ztcol))[-(1:length(fullzlim))]
  at = match(ztat,ztcol)
  plotColorLegend(x0,x1,y0,y1,col,at=at,legend=leg,title=title)
}

plotColorLegend = function(x0,x1,y0,y1,col,at,legend,title=NULL){
	xpd = par(xpd=TRUE)
	y = seq(grconvertY(y0,'nfc','user'),grconvertY(y1,'nfc','user'),length.out = length(col)+1)
	rect(grconvertX(x0,'nfc','user'),y[-length(y)],grconvertX(x0+(x1-x0)*0.25,'nfc','user'),y[-1],col=col,border = NA)
	at = y[at]+(y[2]-y[1])/2
	text(grconvertX(x0+(x1-x0)*0.3,'nfc','user'),at,legend,adj=c(0,0.5))
	if(!is.null(title)){
	  #text(grconvertX(x1,'nfc','user'),y[length(y)],title,adj=c(1,-0.5))
	  text(grconvertX(x0,'nfc','user'),y[length(y)],title,adj=c(0,-0.5))
	}
	par(xpd=xpd)
}

mergePNG2PFD = function(dir=NULL,fls=NULL,pdfout,...){
  require(png)
  pdf(pdfout,...)
  
  if(is.null(fls))
    fls = sort(list.files(dir,full.names = T,pattern = "*.png"))
  
  for(f in fls) {
    print(f)
    pngRaster = readPNG(f)
    plot.new()
    par(xpd=NA)
    rasterImage(pngRaster, grconvertX(0,'ndc','user'), grconvertY(0,'ndc','user'),grconvertX(1,'ndc','user'),grconvertY(1,'ndc','user'))
  }
  dev.off()
}



write.fasta = function(seq,snames = names(seq),file = NULL,len=60){
	if(!is.null(file)) sink(file)
	for(i in 1:length(seq)){
		j = 1
		cat('>',snames[i],'\n',sep='')
		repeat{
			cat(substr(seq[i],j,j+len-1),'\n',sep='')
			j = j + len
			if(j > nchar(seq[i])) break
		}
	}
	if(!is.null(file)) sink(NULL)
}

barplotWithText = function(x,t=x,adj=c(0.5,1.1),srt=0,text.y=x,ylim=c(0,max(x)),...){
	if(length(ylim)==1)
		ylim = c(0,max(x)*ylim)
	b=barplot(x,ylim=ylim,...)
	text(b,text.y,t,adj=adj,srt=srt)
}
#barplotWithText(1:4,adj=c(0.5,-0.1),ylim=c(0,4.2))

number2bin = function(v,n){
	o = order(v)
	j = 1
	for(i in 1:length(v)){
		if(j < i/length(v)*n)
			j = j + 1
		v[o[i]] = j
	}
	v
}

compNt = function(s1,s2){
	na = which(!(is.na(s1) | is.na(s2)) & nchar(s1) == nchar(s2))
	nt = c('A','T','G','C')
	s1 = lapply(strsplit(toupper(s1),''),function(x){x[!(x %in% nt)] = NA;x})
	s2 = lapply(strsplit(toupper(s2),''),function(x){x[!(x %in% nt)] = NA;x})
	r = rep(NA,length(s1))
	r[na] = sapply(na,function(i)all(s1[[i]] == s2[[i]],na.rm=T))
	r
}

addIsCogingByEnsGTF = function(ens.gtf,as){
	t=Sys.time()
	
	cds = read.table(ens.gtf,sep = '\t',header=FALSE,comment.char = '#')[,c(1,3:5,7)]
	cds = cds[cds[,2] == 'CDS',-2]
	colnames(cds) = c('chr_id','start','stop','strand')
	
	seqinf = Seqinfo(unique(c(cds$chr_id,as$seg$chr_id)))
	cds = GRanges(cds$chr_id ,IRanges(cds$start,cds$stop),cds$strand,seqinfo = seqinf)
	cds = reduce(cds)
	
	seg.r = GRanges(as$seg$chr_id,IRanges(as$seg$start,as$seg$stop),ifelse(is.na(as$seg$strand),'*',ifelse(as$seg$strand==1,"+","-")),seqinfo = seqinf)
	overlap = findOverlaps(seg.r,cds,type='any',ignore.strand = FALSE,select='all')@from
	cat('map1: ',as.numeric(Sys.time()-t, units = "secs"),"\n")
	include = findOverlaps(seg.r,cds,type="within",ignore.strand = FALSE,select='all')@from
	cat('map2: ',as.numeric(Sys.time()-t, units = "secs"),"\n")
	as$seg$cod = 'n'
	as$seg$cod[overlap] = 'p'
	as$seg$cod[include] = 'c'
	cod.genes = as$seg$gene_id[as$seg$cod != 'n']
	as$seg$cod.gene = as$seg$gene_id %in% cod.genes
	as
}

calcMeanCols = function(d,f,FUN=base::mean,verbose=FALSE){
	u = sort(unique(as.character(f)))
	r = matrix(ncol=length(u),nrow=nrow(d))
	colnames(r) = u
	rownames(r) = rownames(d)
	for(j in 1:length(u)){
		i = u[j]
		if(verbose) cat('\r',j,' from ',length(u),'; ',i,'      ')
		r[,i] = apply(d[,f==i,drop=F],1,FUN,na.rm=TRUE)
	}
	r
}

revComplement = function(s,rev=TRUE){
  a = c(a='t',g='c',r='y',s='s',w='w',k='m',b='v',d='h','-'='-',n='n')
  a = c(a,setNames(names(a),a))
  a = c(a,setNames(toupper(a),toupper(names(a))),c('['=']','/'='/',']'='[','>'='>'))
  if(rev)
    r = sapply(strsplit(s,''),function(x)paste(a[rev(x)],collapse = ''))
  else
    r = sapply(strsplit(s,''),function(x)paste(a[    x ],collapse = ''))
  r[is.na(s)] = NA
  r
}

plotPanelLetter = function(l,cex=1.2,adj=c(0,1.1),...){
	x=grconvertX(0,from='nfc',to='user')
	y=grconvertY(1,from='nfc',to='user')
	text(x=x,y=y,labels=l,adj=adj,font=2,cex=cex,xpd=NA)
}

function(d,t=NULL,digits=2,text.col=NULL,xaxlab=rownames(d),yaxlab=colnames(d),...){
	if(is.null(t))
		t = round(d,digits = digits)
	x = 1:nrow(d)
	y = 1:ncol(d)
	image(x,y,d,xaxt='n',yaxt='n',...)
	if(is.null(text.col))
		text.col = 'black'
	text(rep(x,times=length(y)),rep(y,each=length(x)),t,col=text.col)
	if(!is.null(xaxlab))
		axis(1,x,xaxlab)
	if(!is.null(yaxlab))
		axis(2,y,yaxlab)
}


plotArea = function(x,p,col,sd.mult=2,new=FALSE,ylim=NULL,xlim=range(x),area.transp=0.2,type='l',area.den=-1,...){
	#p should contain either mean and sd
	#or mean, lower and upper bounds
	o = order(x)
	x = x[o]
	p = p[o,,drop=F]
	na = !is.na(p[,1])
	x = x[na]
	p = p[na,,drop=F]
	if(ncol(p)==2)
		yp = c(p[,1]-p[,2]*sd.mult,rev(p[,1]+p[,2]*sd.mult))
	else
		yp = c(p[,2],rev(p[,3]))
	if(new){
		if(is.null(ylim))
			ylim = range(yp,na.rm=T)
		plot(1,t='n',xlim=xlim,ylim=ylim,...)
	}
	col.pol = col2rgb(col,alpha=TRUE)/255
	col.pol = rgb(col.pol[1],col.pol[2],col.pol[3],col.pol[4]*area.transp) 
	polygon(c(x,rev(x)),yp,col=col.pol,border=NA,den=area.den)
	lines(x,p[,1],col=col,type=type,...)
}

reorderClustersBySize = function(clusts){
  stop('use renameClustsBySize instead')
  t = sort(table(clusts),decreasing=T)
  clusts. = clusts
  for(c in 1:length(t))
    clusts[clusts.==names(t)[c]] = c
  clusts
}


plotLine = function(x,y=NULL,cor.method='pearson',line.col='red',leg.pos='topright',line.lwd=1,plot.ci=FALSE,ci.transparency=0.3,line.on.top=TRUE,new=TRUE,...){
	if(is.null(y)){
		y = x[,2]
		x = x[,1]
	}
  if(new)
    plot(x,y,t='n',...)
  if(line.on.top)
    points(x,y,...)
  f = !is.na(x) & !is.na(y) & !is.infinite(x) &  !is.infinite(y)
	x = x[f]
	y = y[f]
	o = order(x)
	y = y[o]
	x = x[o]
	x = as.numeric(x)
	y = as.numeric(y)
	m = lm(y~x)
	ci=predict.lm(m,interval='confidence')
	if(plot.ci){
		c=col2rgb(line.col)
		polygon(c(x,rev(x)),c(ci[,2],rev(ci[,3])),border=NA,col=rgb(c[1],c[2],c[3],ci.transparency*255,maxColorValue = 255))
	}
	lines(x,ci[,1],col=line.col,lwd=line.lwd)
	if(!line.on.top)
		points(x,y,...)
	if(tolower(cor.method)==substr('spearman',1,nchar(cor.method))){
		x = rank(x)
		y = rank(y)
		cor.method = 'pearson'
	}
	c = cor.test(x,y,m=cor.method)
	ci = round(c$conf.int,2)
	leg=paste(toupper(substr(cor.method,1,1)),'CC=',round(c$estimate,2),' [',ci[1],',',ci[2],']\npv=',format(c$p.value,digits=2,scientific=TRUE),'\nN=',length(x),sep='')
	text(grconvertX(0.01,'npc','user'),grconvertY(0.99,'npc','user'),leg,col=line.col,adj=c(0,1),font=2)
}

plotClustBoxes=function(v,cl,...){
	t = table(cl)
	for(i in 1:length(t)){
		boxplot(v[cl==i,],las=2,main=paste('c',i,' (',t[i],')',sep=''),...)
	}
}

renameClustsBySize = function(x){
	t = table(x)
	n = names(t)
	o = order(t,decreasing=T)
	r = x
	for(i in 1:length(o))
		r[x == n[o[i]]] = i
	r
}

normRows=function(d){
	d = sweep(d,1,apply(d,1,mean,na.rm=T),'-')
	sweep(d,1,apply(d,1,sd,na.rm=T),'/')
}

getPal=function(col=c('blue','white','red'),n=200){
	r = character(n)
	if(n == 1){
		i = floor(length(col)/2)
		return(.getPal(col[i],col[i+1],0.2))
	}
	for(i in 0:(n-1)){
		cinx = i/(n-1)*(length(col)-1)+1
		if(cinx == floor(cinx))
			r[i+1] = col[cinx]
		else{
			rate = cinx-floor(cinx)
			cinx = floor(cinx)
			r[i+1] = .getPal(col[cinx],col[cinx+1],1-rate)
		}
	}
	r
}

.getPal = function(c1,c2,r){
	c1 = t(col2rgb(c1,alpha=TRUE)/255)
	c2 = t(col2rgb(c2,alpha=TRUE)/255)
	return(rgb((c1*r+c2*(1-r)),alpha=r*c1[4]+(1-r)*c2[4]))
}

plotLogFCDistr = function(d,thr,bw=0.1,...){
	mn = min(d[,3])
	mx = max(d[,3])
	us = density(d[d[,2]> 0.05,3],bw=bw,from=mn,to=mx)
	sg = density(d[d[,2]<=0.05,3],bw=bw,from=mn,to=mx)
	plot(us,ylim=c(0,max(us$y,sg$y)),xlab='log2FC',lwd=2,col='gray',...)
	lines(sg,lwd=2,col='red',...)
	abline(v=c(-thr,thr),lty=2,col='black')
	legend('topleft',col=c('gray','red','red'),lwd=2,legend=paste(c('not sign','sign',paste('sign & log2FC>',thr,sep='')),' (',c(sum(d[,2]>0.05),sum(d[,2]<=0.05),sum(d[,2]<=0.05 & abs(d[,3])>=thr)),')',sep=''))
}

plotHM = function(data,main='',method='sp',norm=T,col=c('blue','white','red'),col.count=200,reorder=1:ncol(data),hclust.method="complete",...){
	if(norm)
		data = normRows(data)
	c = cor(data,use='pair',method=method)
	cols = getPal(col,col.count)
	if(norm){
		r = round(range(c)*100)+100
		cols=cols[r[1]:r[2]]
	}
	heatmap(c,symm=T,col=cols,main=paste(main,' (',dim(data)[1],')',sep=''),Rowv=reorder,reorderfun=function(d,w){reorder(d,w,agglo.FUN=mean)},
					distfun=function(d){as.dist(1-d)},hclustfun=function(d){hclust(d,method=hclust.method)},...)
	cx = grconvertX(c(0.03,0.13), from='ndc',to='u')
	cy = grconvertY(c(0.95,0.75), from='ndc',to='u')
	lg = round((0:3)/3*(max(c)-min(c))+min(c),2)
	lg[length(lg)] = '1.00'
	plotrix::color.legend(cx[1],cy[1],cx[2],cy[2],legend=lg,rect.col=cols,gradient="y",align='rb')
}

calcPairCor = function(x,y,...){
	r = numeric(nrow(x))
	for(i in 1:nrow(x))
		r[i] = cor(x[i,],y[i,],...)
	r
}
scaleTo = function(x,from=0,to=1,minx=min(x,na.rm=TRUE),maxx=max(x,na.rm=TRUE),fraction=1){
	x = (x-minx)/(maxx-minx)
	x*(to-from)*fraction + from + (to-from)*(1-fraction)/2
}

getReadCoverage = function(bams,chr,start,end,strand=NA){
	require(GenomicAlignments)
	if(start>end){
		t = start
		start = end
		end=t
	}
	r = list(x = start:end,cov = NULL,juncs=NULL)
	param = ScanBamParam(flag=scanBamFlag(isMinusStrand=strand==-1),which=GRanges(chr, IRanges(start, end)))
	i = 1
	for(b in bams){
		cat('\r',i,'     ')
		i = i + 1
		bam = readGAlignments(b,param = param)
		cov=coverage(bam)[[chr]][start:end]
		juncs = as.data.frame(summarizeJunctions(bam))
		rownames(juncs)=paste(juncs$seqnames,juncs$start,juncs$end,sep='-')
		if(is.null(r$cov)){
			r$cov=cov
			r$juncs=juncs
		}else{
			r$cov = r$cov + cov
			cmn = intersect(rownames(juncs),rownames(r$juncs))
			r$juncs[cmn,'score'] = r$juncs[cmn,'score'] + juncs[cmn,'score']
			r$juncs = rbind(r$juncs,juncs[setdiff(rownames(juncs),rownames(r$juncs)),])
		}
	}
	invisible(r)
}


#' Title
#'
#' @param cov list with three elements: x (coordinates), cov (coverage) and juncs (start,end,score)
#'
#' @return
#' @export
#'
#' @examples
plotReadCov = function(r,min.junc.cov=0,plot.junc.only.within=FALSE,ylim=NULL,xlim=range(r$x),reverse=FALSE,junc.col='red',junc.lwd=3,...){
	f = r$x >= xlim[1] & r$x <=xlim[2]
	r$x = r$x[f]
	r$cov = r$cov[f]
	r$juncs = r$juncs[r$juncs$start <= xlim[2] | r$juncs$end >=xlim[1],]
	
	start = r$x[1]
	end = r$x[length(r$x)]
	r$cov[c(1,length(r$cov))] = 0
	if(is.null(ylim))
		ylim = c(0,max(r$cov,r$juncs$score))
	if(reverse)
		xlim=rev(xlim)
	plot(r$x,r$cov,t='n',ylim=ylim,xlim=xlim,...)
	polygon(r$x,r$cov,col = 'gray',border=NA)
	if(nrow(r$juncs)>0)
		for(i in 1:nrow(r$juncs))
			if(r$juncs$score[i] >= min.junc.cov & (!plot.junc.only.within || (r$juncs$start[i] > min(xlim) & r$juncs$end[i]< max(xlim))))
				plotArc(r$juncs$start[i],r$juncs$end[i],r$juncs$score[i],col=junc.col,lwd=junc.lwd)
}

plotArc = function(from,to,top,n=100,y.base=0,...){
	len = to - from
	x = seq(from=0,to=len,length.out = n)
	y = x*4*top/len - x^2*(4*top/len^2)
	lines(x+from,y+y.base,...)
}

