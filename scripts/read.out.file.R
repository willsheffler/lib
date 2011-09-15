
read.outfile = function(filename, clean=F) {
	tmpfile = p('/tmp/tmp.',runif(1),'.out')
	system(p("grep SCORE ",filename," > ",tmpfile))
	tmp = read.table(tmpfile,header=T,comment.char="#") # this is default comment.char, but make it explicit
	system(p("rm ",tmpfile))
	return( as.data.frame(tmp[,-1] ))
}

read.mtyka.out = function( filename ) {
	head = sapply(read.table(filename,comment.char='',nrows=1,sep=':',strip.white=T),as.character)
	head = c( strsplit(head[1],"[ ]+")[[1]] [-1] , head[-1] )
	body = read.table(filename,comment.char='',skip=1,header=F)[-19]
	names(body) = head
	return(as.data.frame(body))
}

read.boinc.sc = function( filename , keepcols = NULL, skip=1) {
	head = sapply(read.table(filename,comment.char='',skip=skip,nrows=1,sep=':',strip.white=T),as.character)
	head = c( strsplit(head[[2]],"[ 	]+")[[1]], 'userid' )
	body = read.table(filename,comment.char='',skip=2,header=F)
	body = body[,-1]
	names(body) = head
	if( !is.null(keepcols) )
		body = body[,keepcols]
	return( body )
}

read.fasc = function( filename , keepcols = NULL, skip=0) {
	head = sapply(read.table(filename,comment.char='',skip=skip,nrows=1,sep=':',strip.white=T),as.character)
	head = c( strsplit(head[[1]],"[ 	]+")[[1]] )
	body = read.table(filename,comment.char='',skip=skip+1,header=F)
	names(body) = head
	if( !is.null(keepcols) )
		body = body[,keepcols]
	return( body )
}

read.fasc.files = function(directory='.',pattern=NULL) {
	if(is.null(pattern))
		pattern = '.*[.]fasc'
	scores = list()
	for( file in dir(directory,pattern=pattern) ) {
		pp('reading ',file)
		scores[[file]] = read.fasc(p(directory,'/',file))		
	}
	return( scores )
}

