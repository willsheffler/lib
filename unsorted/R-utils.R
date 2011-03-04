#################################################################
#################################################################

#source("~/.env.R")
#source(paste(SCRIPTSDIR,'read.out.file.R',sep=''))


#src = function() {
#	source(paste(SCRIPTSDIR,'R-utils.R',sep=''))
#}

##### start analysis

LIBLOC=NULL
library(e1071,lib.loc=LIBLOC)
# decided not to use these because I don't quite get what ROCR is doing
library(gtools,lib.loc=LIBLOC)
library(gdata, lib.loc=LIBLOC)
library(gplots,lib.loc=LIBLOC)
library(ROCR,  lib.loc=LIBLOC)
library(e1071, lib.loc=LIBLOC)

atypes = c("CNH2","COO","CH1","CH2","CH3","aroC","Ntrp","Nhis","NH2O","Nlys","Narg","Npro",
           "OH","ONH2","OOC","S","Nbb","CAbb","CObb","OCbb","Phos","Hpol","Hapo","Haro","HNbb")

resnames <- c("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
              "MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")

RCOLORSDARK = colors()
for(adjective in c('light','pale','grey','gray','snow','white',
				'khaki','peach','blush','lavender','pink','rose','papaya',
				'aliceblue','azure','blanch','burlywood','cornsilk','lemon',
				'honeydew','ivory','cream','lace','seashell','wheat','beige','bisque',
				'linen','gainsboro','moccasin')) {
	RCOLORSDARK = RCOLORSDARK[-grep(adjective,RCOLORSDARK)]
}

#RCOLORSDARK = strsplit("aquamarine4 azure4 bisque4 black blue4 brown4 darkblue darkgreek"," ")

getpr = function(arad) {
	rad = sqrt( arad**2+c( 0.0000000,0.6167577,0.8835883,1.0960643,
	           1.2816093,1.4507516,1.6089340,1.7595542,
				  1.9051524,2.0480757,2.1910429,2.3378414,
	 			  2.4943877,2.6704755,2.8827284,3.1596130,
				  3.5499700,4.1377324,5.0680891,6.5954530 )**2 ) - arad
	return(c(rad[1:16],2.7,3.2,3.7,4.2))
}

ATYPE_LABEL = list()
ATYPE_LABEL[['001']] = "CNH2"
ATYPE_LABEL[['002']] = "COO"
ATYPE_LABEL[['003']] = "CH1"
ATYPE_LABEL[['004']] = "CH2"
ATYPE_LABEL[['005']] = "CH3"
ATYPE_LABEL[['006']] = "aroC"
ATYPE_LABEL[['007']] = "Ntrp"
ATYPE_LABEL[['008']] = "Nhis"
ATYPE_LABEL[['009']] = "NH2O"
ATYPE_LABEL[['010']] = "Nlys"
ATYPE_LABEL[['011']] = "Narg"
ATYPE_LABEL[['012']] = "Npro"
ATYPE_LABEL[['013']] = "OH"
ATYPE_LABEL[['014']] = "ONH2"
ATYPE_LABEL[['015']] = "OOC"
ATYPE_LABEL[['016']] = "NONE"
ATYPE_LABEL[['017']] = "S"
ATYPE_LABEL[['018']] = "Nbb"
ATYPE_LABEL[['019']] = "CAbb"
ATYPE_LABEL[['020']] = "CObb"
ATYPE_LABEL[['021']] = "OCbb"

pad = function(x) {
	x = as.character(x)
	if(nchar(x)==1) return(paste("00",x,sep=""))
	if(nchar(x)==2) return(paste("0" ,x,sep=""))
	return(x)
}


randrange = function(start,end,N) {
	if(N==0) return(c())
	if(end-start<=N) return(start:end)
	return( (order(runif(end-start+1))+start-1)[1:N] )
}

getmodel = function(s) {
	cbind(t(t(s$SV) %*% s$coefs),s$rho)
}

wpred = function(w,X) {
	dim(w) = c(dim(X)[2]+1,1)
   r = as.matrix(X) %*% w[1:dim(X)[2]] - w[dim(X)[2]+1]
	names(r) = rownames(X)
	return(r)
}




res2num  <- 1:20
names(res2num) <- resnames

len <- length

minmax = function(x,integer=F) {
	if( integer )
 		return( c(floor(min(x)),ceiling(max(x))) )
	else
 		return( c(min(x),max(x)) )
}

p = function(...) paste(...,sep='')
pp = function(...) print(paste(...,sep=''))


load.compressed = function(file, recompress=T, envir=.GlobalEnv ) { # parent.frame()?
	pp("USE SAVE/LOAD compress=T!!!!!!!!!!!!!")
}

#endswith = function(target, pattern) {
#	
#}

basename = function( str ) {
	sub('^.*[/]','',str)
}

pathname = function( str ) {
	sub('[^/]*$','',str)
}


frac.same.consecutive = function(x) {
	sum( x[-1]==x[-length(x)] ) / len(x)
}
frac.diff.consecutive = function(x) {
	sum( x[-1]!=x[-length(x)] ) / ( len(x)-1 )
}

clean.up.data.frame = function(frame,namesfromcol=1,sort=T,prune.loners=T) {

	rn = frame[,namesfromcol]
	rn = tapply(1:len(rn),rn,min)
	if(sort) {
		rn = rn[ order( names(rn) ) ] # do this first to ensure .. bla bla
		#rn = rn[ order( gsub('/(.*?/)','',names(rn)) ) ] # quicksort would be best
	}	#																  # for this interdigitated list
	frame = frame[rn,]
	rownames(frame) = names(rn)
	return(frame)
	
}

isectlist = function( listofstuff, func=I ) {
	isection = func(listofstuff[[1]])
	for( ii in 2:len(listofstuff)) {
		isection = intersect(isection,func(listofstuff[[ii]]))
	}
	return(isection)	
}

homogenize.data.frames = function(framelist,...) {
	for( n in names(framelist) )
		framelist[[n]] = clean.up.data.frame( framelist[[n]] )
	sharednames = isectlist( framelist, rownames )
	return( sharednames )
}

read.table.compressed = function(file,recompress=T,...) {
  n = nchar(file)
  was.compressed = F
  if( substr(file,n-2,n) == ".gz" ) {
    p("uncompressing ",file)
    system(p("gzip -d ",file))
    file = substr(file,1,n-3)
    was.compressed = T
  }
  p('reading ',file)
  tmp = read.table(file, ... )
  if( recompress ){
    p('recompressing ',file)
    system( p("gzip ",file) )
  }
  return( tmp )
}


outerpairs = function(x,y) {
  M = max(c(x,y))+1
  p = c()
  for(i in x) {
    for(j in y) {
      if( i < j ) {
        p = append(p,M*i+j)
      } else {
        p = append(p,M*j+i)
      }
    }
  }
  p = unique(p)
  print(p)
  l = list()
  for(e in p) {
    l[[ len(l)+1 ]] = c( floor(e/M), e%%M )
  }
  return( l )
}

maxima = function(x,thresh=1e-6,notails=T) {
  m <- c()
  if( x[1] > x[2] )
    m <- append(m,1)
  for(ii in 2:(len(x)-1) ) {
    #if(x[ii] >  thresh+x[ii-1] & x[ii] >=        x[ii+1]  |
    #   x[ii] >=        x[ii-1] & x[ii] >  thresh+x[ii+1])
    if(x[ii] >  thresh+x[ii-1] & x[ii] >  thresh+x[ii+1])
      m <- append(m,ii)
  }
  if( x[length(x)] > x[length(x)-1] )
    m <- append(m,length(x))
  m = unique(c(m,which(x==max(x))[1]))
  return(m)
}

minima = function(x,thresh=1e-6,notails=F) {
  return( maxima(-x,notails) )
}



find.limits.single = function(data,sdth=5,yrange=NULL,notails=T) {  
  #print("find.limits.single")
  ub <- data[maxima(data,notails=notails)]
  lb <- data[minima(data,notails=notails)]
  m  <- median(data)
  sd <- mad(data)
  ub <- ub[ abs(ub-m) <= sdth*sd ]
  lb <- lb[ abs(lb-m) <= sdth*sd ]
  if(length(ub)==0) ub <- NA
  if(length(lb)==0) lb <- NA
  ub <- max(ub)
  lb <- min(lb)
  if( !is.na(ub) & !is.na(lb) ) {
    if( is.null(yrange) )
      return( c(lb,ub) )
    else if( ub-lb > yrange )
      return( c(lb,ub) )
    else
      return( mean(c(lb,ub))+ yrange/2*c(-1,1) )
  } else if( is.na( lb )  ) {
    if( is.null(yrange) ) 
      if( ub < 0 )
        return( c(ub,-ub/4) )
      else
        return( c(ub-1,ub ) )
    return( c( ub-yrange*19/20,ub+yrange/20 ) )
  } else if( is.na( ub )  ) {
    if( is.null(yrange) )
      if( lb < 0 )
        return( (c(lb,-lb/4)) )
      else
        return( c(lb,lb+1))
    return( c( lb-yrange/20,lb+yrange*19/20 ) )
  } else {
    return( median(data) + yrange/2*c(-1,1) )
  }
}


numericintegral <- function(x,y) {
  dx = x[2:len(x)]-x[1:(len(x)-1)]
  #dy = x[2:len(y)]-x[1:(len(y)-1)]
  dx = c(dx[1],dx)
  #dy = c(dy[1],dy)
  ipos = x>0
  ineg = x<0
  #print(paste(len(dx),len(y)))  
  sumpos = sum(dx[ipos]*y[ipos])
  sumneg = sum(dx[ineg]*y[ineg])
  sum = sum(dx*y)
  return(list(sum=sum,sumpos=sumpos,sumneg=sumneg))
}


srcutil = function() source('~/scripts/R-utils.R')


mysource <- function(file) {
  if( substr(file,1,6) == "/Users" )
    file <- paste("/users",substr(file,7,nchar(file)),sep='')
  
  source(file)
}
myload <- function(file) {
  if( substr(file,1,6) == "/Users" )
    file <- paste("/users",substr(file,7,nchar(file)),sep='')
  load(file)
}

numdiff = function(x,y){
  dx = c( x[2]-x[1] , x[2:len(x)] - x[1:(len(x)-1)] )
  dy = c( y[2]-y[1] , y[2:len(y)] - y[1:(len(y)-1)] )
  return ( dy/dx )
}

plotden <- function(x, y, z=NULL, cond=function(x) 1:length(x), title="", bw=1,
         rocfrac=1, nroc=2, adjust=0.1, xlim=c(0,1), ... ) {  
  if(!is.null(z)) {
    print('asldfkjhasldkfjlksd')                                        #r <- calc.roc(x,z,n=50)
    dx <- density(x,adjust=adjust)
    dy <- density(y,bw=dx$bw)    
    dz <- density(z,bw=dx$bw)
    ymax <- max(dx$y,dy$y,dz$y)
    plot(dx,ylim=c(0,ymax),main=paste(title),xlim=xlim,...)
    lines(dy,col=2,...)
    lines(dz,col=3,...)
  } else {
    #if( min(x,y) < xlim[1] ) # allow expansion of xlim but not contraction
    #  xlim[1] = min(x,y)
    #if( max(x,y) > xlim[2] )
    #  xlim[2] = max(x,y)
    r = NA
    try( { 
	r    <- getroc(x,y) })
    dx   <- density(x,adjust=adjust)
    dy   <- density(y,bw=dx$bw)    
    ymax <- max(dx$y,dy$y)    
    xmx  <- quantile(c(dx$x,dy$x),0.5)
    plot(0,type='n',ylim=c(0,ymax),xlim=xlim,main=paste(title,round(r,4)),...)
    lines(dx,col=1,...)    
    lines(dy,col=2,...)    
    return(r)
  }  
}

getroc <- function(x,y) {
  predictions = c(x,y)
  labels      = c(rep(T,length(x)),rep(F,length(y)))
  pred <- prediction(predictions, labels)
  perf <- performance(pred, measure = "auc")
  return(slot(perf,"y.values")[[1]])
}



first <- function() {
  if( !(2 %in% dev.list()) )
    dev.x11()
  if( !(2 %in% dev.list()) ) {
    print("can't switch to dev 2!")
    return
  }
  dev.set(2)
}
second <- function() {
  if( !(3 %in% dev.list()) )
    dev.x11()
  if( !(3 %in% dev.list()) ) {
    print("can't switch to dev 3!")
    return
  }
  dev.set(3)
}
third <- function() {
  if( !(4 %in% dev.list()) )
    dev.x11()
  if( !(4 %in% dev.list()) ) {
    print("can't switch to dev 4!")
    return
  }
  dev.set(4)
}


calc.joint.bins <- function(bp,fieldlist,namesubs,fsep=".",nsep="") {
  if(missing(namesubs))
    sub <- list()
  bins <- list()
  bins[['nobins']] <- as.factor(rep('nobins',dim(bp)[1]))
  for(ii in 1:length(fieldlist)) {
    fields <- fieldlist[[ii]]
    name <- paste(fieldlist[[ii]],collapse=nsep)
    for(sub in namesubs)
      name <- sub(sub[1],sub[2],name)
    print(paste(c("calculating joint factor on fields:",fields),collapse=" "))
    jointfield <- as.character(bp[[fields[1]]])
    if(length(fields)>1) {
      for(jj in 2:length(fields)) {
        field <- fields[jj]
        jointfield <- paste(jointfield,as.character(bp[[field]]),sep=fsep)
      }
    }
    bins[[name]] <- as.factor(jointfield)
  }
  return(bins)  
}

calc.frac.unsatisfied <- function(unsatisfied,bins) {
  uf <- list()
  for(bin in names(bins)) {
    ct <- tapply(unsatisfied, bins[[bin]], length)
    uc <- tapply(unsatisfied, bins[[bin]], sum)
    uf[[bin]] <- uc/ct
  }
  return(uf)
}

calc.pred <- function(bins,predictor,satisfied=F,rm.satisfied=F) {
  pred <- list()
  for(bin in names(bins)) {
    tmp <- bins[[bin]]
    levels(tmp) <- predictor[[bin]]
    tmp <- as.numeric(as.character(tmp))
    if(rm.satisfied)
      tmp <- tmp[!satisfied]
    else if(any(satisfied))
      tmp[satisfied] <- 1
    pred[[bin]] <- tmp
  }
  return(pred)
}


# function to plot roc curves
# labels must be true and false
calc.roc <- function(posval,negval,n=10,highgood=T) {
  pred  <- c(posval,negval)
  label <- c(rep(T,length(posval)),rep(F,length(negval)))
  tpr <- 0:(n-1)/(n-1)
  fpr <- 0:(n-1)/(n-1)
  q <- quantile(pred,0:(n-1)/(n-1))  
  for(ii in 2:(n-1)) {
    x <- sum(pred<=q[ii] & !label)/sum( !label )
    y <- sum(pred<=q[ii] &  label)/sum(  label )
    fpr[ii] <- x
    tpr[ii] <- y
  }
  if(highgood) {
    tpr <- rev(1-tpr)
    fpr <- rev(1-fpr)
  }
  score <- 0
  for(ii in 2:n) {
    #print(paste(ii,score, mean(c(tpr[ii-1],tpr[ii])),abs(fpr[ii]-fpr[ii-1])))
    score <- score + mean(c(tpr[ii-1],tpr[ii]))*abs(fpr[ii]-fpr[ii-1])
  }
  return(list(y=tpr,x=fpr,score=score))
}

calc.roc.categories <- function(posvals,negvals,n=10,highgood=T) {
  if(any(names(posvals)!=names(negvals)))
    print("ERROR calc.roc.categories names for pos and neg don't match")
  roc <- list()  
  for(bin in names(posvals)) {
    roc[[bin]] <- calc.roc(posvals[[bin]],negvals[[bin]],n,highgood)
  }
  return(roc)
}



plot.roc <- function(roc,title="Some ROC Curves",colors=NULL,file=NULL,legx=0.50,legy=0.35) {
  rocscores <- format(round(as.numeric(lapply(roc,function(x) x$score)),3))
  if(is.null(colors))
    colors = rainbow(length(roc))
  plot (c(0,1),c(0,1) , col=1 , type='n',main=title,
        xlab="Fraction False Positives",ylab="Fraction True Positives")
  for(ii in 1:length(roc)) 
    lines(roc[[ii]], col=colors[ii])
  legend(legx,legy,c(paste(rocscores,'on:',names(roc))), col=c(colors),
         lwd=c(rep(3,length(roc))))
  if(!is.null(file))
    dev.print(file=file,device=postscript,horizontal=F)       
}



mydensity <- function(x) {
  if(length(x)==1) return(list(x=c(x,x),y=c(0,1)))
  else             return(density(x))
}

gt  <- function(a,b) return( a >  b )
gte <- function(a,b) return( a >= b )
lt  <- function(a,b) return( a <  b )
lte <- function(a,b) return( a <= b )

padzeros <- function(x,n=2) {
  x <- as.character(x)
  y <- x
  for(ii in 1:length(x))
    y[ii] <- paste( c(rep("0",(n-nchar(x[ii]))) , x[ii]) , sep="" , collapse='' )
  return(y)
}

calcbins <- function(data,breaks){
  bin <- rep(NA,length(data))
  for(ii in 1:(length(breaks)-1)) {
    lcmp <- lte
    rcmp <- lt
    if( ii == 1) {
      lcmp <- lte
    } else if( breaks[ii] == breaks[ii+1] ) {
      lcmp <- lte
    } else if( breaks[ii-1] == breaks[ii] ) {
      lcmp <- lt
    }
    if( ii == length(breaks)-1 ) {
      rcmp <- lte
    } else if(breaks[ii]==breaks[ii+1]) {
      rcmp <- lte
    } else if(breaks[ii+1]==breaks[ii+2]) {
      rcmp <- lt
    }
    #lbl <- paste(format(breaks[ii],  digits=3),'-',
    #             format(breaks[ii+1],digits=3),sep='')
    lbl <- padzeros(ii,2)
    #lbl <- padzeros(ii,nchar(as.character(length(breaks)-1)))
    bin[ lcmp(breaks[ii],data) & rcmp(data,breaks[ii+1]) ] <- lbl
  }
  bin <- as.factor(bin)
  return(bin)
}
#tmp <- calcbins(pdbe$sasafrac,c(0,0,0.5,1,1))

calcbinsbyfactor <- function(data,fac,breaklist) {
  bin <- c()
  for(l in levels(fac)){
    b <- calcbins(data[fac==l],breaklist[[l]])
    bin[fac==l] <- as.character(b)
  }
  bin <- as.factor(bin)
  return(bin)
}

# TODO: write joint bins function
calcjointbins <- function(data,breaks) {
  
}

#calcsasabins <- function(sasafrac) {
#  return(calcbins(sasafrac,c(0,0,0.02,0.05,0.09,0.15,0.22,0.30,0.45,0.60,1.0)))
#}


quantilebreaks <- function(data,nbins) {
  # TODO!!! make this take big modes as special case rather than just 0!
  N <- length(data)
  if(sum(data==0)/N > 1/nbins) {
    return(c(0,0,quantile(data[data!=0],1:(nbins-1)/(nbins-1))))
  } else {
    return(quantile(data,0:nbins/nbins))
  }
}

quantilebreaksbyfactor <- function(data,fac,nbins) {
  breaks <- as.list(1:length(levels(fac)))
  for(ii in 1:length(levels(fac))){
    breaks[[ii]] <- quantilebreaks(data[fac==levels(fac)[ii]],nbins)
    breaks[[ii]][1] <- min(data)
    breaks[[ii]][length(breaks[[ii]])] <- max(data)
  }
  names(breaks) <- levels(fac)
  return(breaks)
}

#jointquantilebreaks <- function(data,nbins) {
#  n <- dim(data)[,2]
#  breaks <- as.list()
#  breaks1 <- quantilebreaks(data,nbins1)
#  breaks2 <- quantilebreaks(data,nbins2)
#  return(list(breaks1,breaks2))
#}

# TODO: write joint breaks function
# should specify desired # per bin and some kind of weight to prefer more
# bins for one "factor" vs. another
# do this based on quantiles of whitened data!

# X <- cbind(pdbe$nb+runif(dim(pdbe)[1])-0.5,pdbe$sasafrac)
# pc <- prcomp(X)
# C <- t(t(X)-pc$center)
# S <- t(t(C)/pc$sdev)
# R <- S %*% pc$rotation

# nb <- pdbe$nb + runif(length(pdbe$nb))
# sf <- pdbe$sasafrac
# nb <- ( nb - mean(nb) ) / sd( nb )
# sf <- ( sf - mean(sf) ) / sd( sf )

# plot( nb , sf , pch='.')

list2array <- function(l) {
  if(is.numeric(l))
    return(l)
  ii <- 1
  while(is.null(l[[ii]]))
    ii <- ii+1        
  d <- c(dim(l),length(l[[ii]]))
  n <- append(dimnames(l),list(names(l[[ii]])))
  a <- array(NA,d)
  dimnames(a) <- n
  if(length(dim(l))==1) {
    for(ii in 1:d[1])
      if(!is.null(l[[ii]]))
        a[ii,] <- l[[ii]]
    return(a)
  } else if(length(dim(l))==2) {
    for(ii in 1:d[1])
      for(jj in 1:d[2])
        if(!is.null(l[[ii,jj]]))
          a[ii,jj,] <- l[[ii,jj]]
    return(a)
  } else if( length(dim(l))==3 ) {
    for(ii in 1:d[1])
      for(jj in 1:d[2])
        for(kk in 1:d[3])
          if(!is.null(l[[ii,jj,kk]]))
            a[ii,jj,kk,] <- l[[ii,jj,kk]]
    return(a)
  } else {
    print(paste('list dim',length(dim(l)),'not handled yet'))
  }
}

write.data.frame <- function(data,file,pretty=F,meta,verbose=F) {
  if(verbose) print(paste("writting data frame",file))
  n <- colnames(data)
  top <- as.data.frame(as.list(n))
  colnames(top) <- n
  for(ii in 1:length(n))
    top[,ii] <- as.character(top[,ii])
  data <- rbind(top,data)
  if(pretty) {
    for(ii in 1:length(n)) {
      #print(paste("padding column ",ii," (slow; find better way!)"))
      #l <- max(nchar(data[,ii]))
      data[,ii] <- format(data[,ii],justify='left') #sapply(data[,ii],function(x) padright(x,l))
    }
  }
  if(missing(meta)) {
    write.table(data,file=file,quote=F,row.names=F,col.names=F,append=F)
  } else {
    write(meta,file=file)
    write.table(data,file=file,quote=F,row.names=F,col.names=F,append=T)
  }
}

nulleplot = function() {
  plot( 0,0, type='n', xlim=c(3,6),ylim=c(-1,1)/3 )
}

plot.with.densities = function(x,y,MARGIN=0.2,ADJUST=0.5,...) {
	XSIZE = max(x)-min(x)
	YSIZE = max(y)-min(y)
	XMARGIN = XSIZE*MARGIN
	YMARGIN = YSIZE*MARGIN
	XLIM = c(min(x)-XMARGIN,max(x))
	YLIM = c(min(y)-YMARGIN,max(y))
	plot(x,y,xlim=XLIM,ylim=YLIM,...)
	dx = density(x,adjust=ADJUST)
	dy = density(y,adjust=ADJUST)
	dx$y = dx$y * YMARGIN*0.8 / max(dx$y)
	dy$y = dy$y * XMARGIN*0.8 / max(dy$y)
	lines(dx$x,YLIM[1]+dx$y,col=4)
	lines(XLIM[1]+dy$y,dy$x,col=4)
}

# x = runif(1000)
# y = rnorm(1000)
# plot.with.densities(x,y)
