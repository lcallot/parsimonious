
# A small function to compute the degrees of freedom for the ridge regresson
.ridge.df<-function(x,lambda){

	x<-(cbind(1,x))
	df.ridge <- NULL
	for(l in lambda){
		df.ridge<-c(df.ridge,sum(diag(x%*% solve(t(x) %*% x + l*diag(ncol(x)))%*% t(x))))
	}
	return(df.ridge)
}