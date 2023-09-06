\name{BaRanking}
\alias{BaRanking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Barrier permeability quantification
}
\description{
Evaluating permeabilities of all barrier segments based on the results of the Barrier Behavior Analysis (\code{BaBA}). Returns an \code{sf MULTILINESTRING} object depicting the permeability index value for each barrier segment with at least \code{min_total_enc} encounters, as well as the number of each type of encounter event that occured near each segment, the total number of encounter events (\code{total_enc}) for each barrier segment, and the number of unique individuals (\code{unique_ind}) that encountered each barrier segment). 
}
\usage{
BaRanking(classification, 
barrier, 
d, 
Barrier_ID,
min_total_enc = 0, 
index_fun = expression(((Bounce + Back_n_forth + Trace + Trapped)/total_enc) * unique_ind),
show_plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classification}{
Classification results from the \code{BaBA} function (obtained with \code{$classification})
}
  \item{barrier}{
  same \code{sf LINESTRING} or \code{MULTILINESTRING} object representing barrier locations used when running the \code{BaBA} function. More see details.
  }
  \item{d}{barrier buffer size used in \code{BaBA} function.}
  \item{Barrier_ID}{The column name from the \code{barrier} object that indicate barrier ID.}
  \item{min_total_enc}{minimum number of encounters to consider when calculating ranking index. More see details.}
  \item{index_fun}{\code{expression} for quantifying permeability of each barrier. More in details.}
  \item{show_plot}{\code{logical}, \code{TRUE} if the barrier map colored by permeability should be printed.}
}
\details{
Classified barrier encounter events will be counted for the barrier segments at which the events occured. If the input barrier contains a few very long segments, users may want to break them into shorter sections so permeabilities will be examined for each smaller section.
 
Depending on animal movement distributions, the barrier segments may only be encountered a few times, or not at all. In this case, there may not be sufficient information to accurately quantify permeabilities of rarely encountered barriers. \code{min_total_enc} is used to filter out barriers with insufficient encounter events from the permeability index calculation.  
 
Depending on management goal and animal ecology/biology, \code{index_fun} can be an expression regarding the total encounter events \code{total_enc}, the number of unique individuals affected \code{unique_ind}, and/or the occurence any of the six barrier behaviors (\code{Bounce, Trace, Back_n_forth, Quick_Cross, Avergae_Movement,} or \code{Trapped}) regarding to each barrier segment. Default is set to be the ratio of counts of altered behavior events (Bounce + Back_n_forth + Trace + Trapped) to counts of total encounter events (\code{total_enc}, weighted by the number of unique individuals that interacted with a given barrier segment. The result of this expression is scaled to range from 0 to 1 using equation: \deqn{function(x){(x-min(x))/(max(x)-min(x))}} \strong{NOTE:} due to the scaling of the permeability index from 0-1, index values are context-specific and cannot be compared across studies or with different input data, other than as a relative indicator of permeability.
}
\value{
an object of class \code{sf} with \code{MULTILINESTRING} geometry containing one row per barrier ID, and columns with information on:
\item{Quick_Cross, Bounce, unknown,
Avergae_movement, Trace, Back_n_forth,
Trapped}{total counts per barrier segment of the corresponding barrier behavior event type.}
\item{unique_ind}{number of unique individuals that encountered the given barrier segment.}
\item{total_enc}{total number of encounters for the given barrier segment.}
\item{index}{permeability ranking index as calculated by \code{index_fun} and transformed to a scale of 0 to 1.}
}
\references{
Xu W, Dejid N, Herrmann V, Sawyer H, Middleton AD. Barrier Behaviour Analysis (BaBA) reveals extensive effects of fencing on wide-ranging ungulates. J Appl Ecol. 2020;00:1–9. https://doi.org/10.1111/1365-2664.13806
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
# load data
data(fences)
data(pronghorn)

# run BaBA function to classify barrier behaviors
results_prong <- BaBA(animal = pronghorn, barrier = fences, d = 110, max_cross = 4)

# rank barriers with default evaluation expression
fence_ranking <- BaRanking(classification = results_prong$classification, barrier = fences, d = 110, Barrier_ID = "FID_Fence0", min_total_enc = 5, show_plot = T)  

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
