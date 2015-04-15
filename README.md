# Graphfilter
Usage: graphfilter [options] [graph.file]

Graphfilter allows you to perform different transformations on the edges
of a graph, like taking the log/-log in any base, rasing to a given power,
taking the exponential, etc, or filtering out some singular or eigenvalues.

Graphs can be read from any table. The only condition is that the nodes
are expected on the first two columns of the table. The location of the
edge values can be specified with option -c.

By extension, it also works with any (non-square) matrix, given that the input abides that condition.
Any matrix element A(i,j) is akin to an edge between nodes `i' and `j'.

Use option -h, --help for help.


Options:<br/>
* General:<br/>
 *  -c n 				Specifies the edges to be at the n-th column<br/>
 *  -h, --help 			Prints this help
 *  -n, --precision n 		Set the number of decimals to use. Does not apply for SVD.
 *  -G, --print-no-diagonal 	When dealing with SV, do no print out the diagonal
 *  -V, --print-svd 		Print matrices U and V from the Singular Value Decomposition of the given graph G: G=USV^t
 *  -q 				Quiet mode: do not print any comment lines.
 *  -Q 				Censor mode: Be quiet and also skip any comment lines from input file. Comment lines are those starting with `#'
 *	-g, --scientific 		Use scientific notation for the output. Does not apply for SVD.
 *	--verbose 			Prints out more detailed info on what's going on.
 *	-X, --non-symmetric 		Do _not_ assume provided matrix is symmetric. By default we assume it is.
* Transformations:
 *	-A, --add w 			Adds w to all edges
 *	-B, --base b  			Specifies the base of logarithms to use
 *	-C, --diag w 			(Re)Sets all diagonal values to w.
 *	-D, --div w 			Divides all edges by w
 *	-E, --exp  			Takes the exp of the edges
 *	-I, --identity  		Do nothing but printing out the edges. Default action.
 *	-L, --log  			Takes the log (base e) of the edges
 *	-N, --log-neg  			Takes the negative log (base e) of the edges
 *	-M, --mul w 			Multiplies all edges by w
 *	-P, --pow n 			Takes the n-power of the edges
 *	-T, --print-triangular 		When dealing with SV, just print out the upper/lower triangular part
 *	-p, --prune-singular-values sv 	Removes all singular values below sv
 *	-v, --singular-values 		Get the singular values of the matrix
 *	-S, --symmetric 		Assume provided matrix is symmetric (Default).

Examples:
*	cat graph | graphfilter -c 11 --log-neg --base 10 > graph-log10

 *	Considers column 11 of file graph as the edge values and takes the -log of them in base 10.

*	cat graph | graphfilter -c 11 --mult 3.5 | graphfilter --exp > graph-exp3.5

 *	The file graph-exp3.5 will contain the exponential of original edges to the power of 3.5, i.e, exp(3.5*edge)
 *	By default graphfilter expects the edges in column 3 (The first column is 1).

* Filter out all singular values below 100 and get the new matrix values:
 *	graphfilter graph2 --prune-singular-values 100

Graphfilter
Copyright (C) Miguel A. Santos, HSC, Toronto, 2009.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )
