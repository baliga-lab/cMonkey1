x <- data.frame()

for ( k in 1:clusterStack$k ) {
	cc <- clusterStack[[ k ]]
	mo <- paste( cc$motif.out$meme.out[[ 1 ]], collapse="<<<<>>>>", sep="<<<<>>>>" )
	xx <- data.frame( k=cc$k, rows=paste(cc$rows,collapse=','), cols=paste(cc$cols,collapse=','), 
		resid=cc$resid, meme_out=mo )
	x <- rbind( x, xx )
}

source( "clusters.R" )