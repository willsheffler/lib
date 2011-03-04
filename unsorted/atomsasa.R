
gm <- function(x)
  exp(mean(log(x)))

ballscoresim <- function(samp) {
  x = rexp(samp*2)/2
  return( 1-x[ x < 1 ][1:samp] )
}
  
  
geommeansim <- function(N,samp=1,generator=runif) {
  x = matrix(generator(N*samp),c(samp,N))
  m = apply(x,1,gm)
  return( m )
}

testnormilization <- function() {
  plot(0,0,type='n',xlim=c(0,1), ylim=c(0,1))
  for( n in 0:8 ) {
    N = 2**n
    dg = density( exp( -rgamma( 100000, N )/N ) )
    de = density( exp( -geommeansim( N, 10000, ballscoresim ) ) )
    
    lines( dg$x, dg$y/max(dg$y) )
    lines( de$x, de$y/max(dg$y) , col=2)
  }
}
