\name{BaRanking}
\alias{BaRanking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Barrier permeability quantification
}
\description{
Evaluating permeabilities of all barrier segments based on the results of the Barrier Behavior Analysis (\code{BaBA}). Returns a dataframe showing the number of each type of encounter events occured near each barrier segment, as well as the total number of encounter events occurred (\code{total_enc}), and the number of unique individuals (\code{unique_ind}) that encountered each barrier segment). 
}
\usage{
BaRanking(classification, 
barrier, 
d, 
min_total_enc = 0, 
index_fun = expression(((Bounce + Back_n_forth + Trace + Trapped)/total_enc) * unique_ind),
show_plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classification}{
Classification results from the \code{BaBA} function (obtained with \code{$classification})
}
  \item{barrier}{
  same \code{SpatialLinesDataFrame} representing barrier locations used when running the \code{BaBA} function. More see details.
  }
  \item{d}{barrier buffer size used in \code{BaBA} function.}
  \item{min_total_enc}{minimum number of encounters to consider when calculating ranking index. More see details.}
  \item{index_fun}{\code{expression} for quantifying permeability of each barrier. More in details.}
  \item{show_plot}{\code{logical}, \code{TRUE} if the barrier maps colored by permeability should be printed.}
}
\details{
Classified barrier encounter events will be counted for the barrier segments at which the events occured. If the input barrier contains a few very long segments, users can break them into shorter sections so permeabilities will be examined for each smaller section.
 
Depending on animal movement distributions, the fence segments may only be encountered for a few times or even are not encountered by animals at all. In this case, we may not have enough information to accurately quantify permeabilities of such barriers. \code{min_total_enc} is used to eliminate these barriers with insufficient encounter events.  
 
Depending on management goal and animal ecology/biology, \code{index_fun} can be an expression regarding the total encounter events \code{total_enc}, the number of unique individuals affected \code{unique_ind}, and the occurence any of the six barrier behaviors (\code{Bounce, Trace, Back_n_forth, Quick_Cross, Avergae_Movement,} or \code{Trapped}) regarding to each barrier segment. Default is set to be the ratio of counts of altered behavior events (Bounce + Back_n_forth + Trace + Trapped) to counts of total encounter events (\code{total_enc}, weighted by numbers of unique individuals interacted with the barriers. Result of this expression is scaled to range from 0 to 1 using equation: \deqn{function(x){(x-min(x))/(max(x)-min(x))}}
}
\value{
an object of class \code{sf} and \code{data.frame} containing one row per fence ID, and columns with information on:
\item{unique_ind}{number of unique individuals that encountered the barrier.}
\item{Quick_Cross, Bounce, unknown,
Avergae_movement, Trace, Back_n_forth,
Trapped}{total counts of the corresponding barrier behavior event type.}
\item{total_enc}{total number of encounters.}
\item{index}{ranking index as calculated by \code{index_fun} and transformed to a scale of 0 to 1.}
}
\references{
Xu, W., Nandintsetsegm D., Sawyer, H., Herrmann, V., Middleton, A. (n.d.). Barrier Behavior Analysis (BaBA) reveals extensive effects of fencing on wide-ranging animals. In Review.
}
\author{
Wenjing Xu \email{wenjing.xu@berkeley.edu} and Valentine Herrmann \email{HerrmannV@si.edu}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(fences)
data(pronghorn)

results_prong <- BaBA(animal = pronghorn, barrier = fences, d = 110, max_cross = 4)

fence_ranking <- BaRanking(classification = results_prong$classification, barrier = fences, d = 110, min_total_enc = 5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
