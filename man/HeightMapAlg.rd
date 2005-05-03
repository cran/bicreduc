\name{HeightMapAlg}
\alias{HeightMapAlg}
\title{HEIGHT MAP ALGORITHM}
\description{
      The HeightMapAlgorithm computes the regions of possible mass support for the
NPMLE for the distribution function of bivariate interval-censored data. 
}
\usage{
HeightMapAlg(R,B)
}
\arguments{
   \item{R}{An n x 4 real matrix of n observation rectangles. Each row corresponds to
       an observation rectangle, represented as (x1,x2,y1,y2). Here (x1,y1) is the
       lower left corner of the rectangle and (x2,y2) is the upper right
       corner of the rectangle. We call (x1,x2) the x-interval and (y1,y2)
       the y-interval of the observation rectangle.
   }
   \item{B}{B describes the boundaries of the observation rectangles (0=open or 1=closed).
       It can be specified in three ways:

       (1) An n x 4 matrix containing 0's and 1's. Each row corresponds to an observation
       rectangle, and is denoted as (cx1, cx2, cy1, cy2). Here cx1 denotes the boundary
       type of x1, cx2 denotes the boundary type of x2, etc. 

       (2) A vector (cx1, cx2, cy1, cy2) containing 0's and 1's. This representation
       can be used if all rectangles have the same type of boundaries. 

	 (3) A vector (c1, c2) containing 0's and 1's. This representation can be
       used if all x- and y-intervals have the same type of boundaries.
       c1 denotes the boundary type of x1 and y1, and c2 denotes the boundary type of 
       x2 and y2.
   }
}
\value{The function returns an m x 4 matrix, containing the maximal intersections,
       that is the areas where the observation rectangles have maximal overlap and where
	the NPMLE can possibly assign mass. Each row (x1,x2,y1,y2) corresponds to a maximal intersection.
}
\references{M.H. Maathuis (2005), "Reduction algorithm for the NPMLE for the distribution
  function of bivariate interval-censored data", Journal of Computational and Graphical
  Statistics (to appear). See also \url{http://www.stat.washington.edu/marloes}}
\author{Marloes Maathuis: \email{marloes@stat.washington.edu}}
\seealso{
}
\examples{
# an example with 6 arbitrarily chosen observation rectangles
R1<-c(1.5, 6.2, 7, Inf)		# first rectangle
R2<-c(2.3, 5, 5.1, Inf)		# second rectangle
R3<-c(3.8, 9.4, 8.3, 10)	# etc...
R4<-c(4.1, Inf, 4.2, 6.7)
R5<-c(7.2, 8.8, 2.7, 9.3)
R6<-c(10, Inf, 1.1, 3.9)
R<-rbind(R1,R2,R3,R4,R5,R6)	# R contains all observation rectangles	
A<-HeightMapAlg(R, c(0,1))	# A contains the maximal intersections
A						
}
\keyword{survival}
\keyword{nonparametric}
\keyword{manip}
}

