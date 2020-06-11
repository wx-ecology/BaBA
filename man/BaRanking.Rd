\name{BaRanking}
\alias{BaRanking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ranking of barriers based on the results of the Barrier Behavior Analysis (\code{BaBA})
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BaRanking(classification = results_pron$classification, 
barrier = fences, 
d = 110, 
min_total_enc = 10, 
index_fun = expression((alt_enc/total_enc)*unique_ind),
show_plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classification}{
Classification results from the \code{BaBA} function (obtained with \code{$classification})
}
  \item{barrier}{
  same \code{SpatialLinesDataFrame} representing barrier locations used when running the \code{BaBA} function.
  }
  \item{d}{barrier buffer size used in \code{BaBA} function.}
  \item{min_total_enc} { minimum number of encounters to consider when calculating ranking index.}
  \item{index_fun}{\code{expression} representing the calculation wished to be used to rank each barrier regarding the number of any encoutner type (\code{Quick_Cross}, \code{Bounce}, \code{unknown}, \code{Average_Movement}, \code{Trace}, \code{Back_n_forth}, \code{Trapped}), the total number of encounters (\code{total_enc} = Bounce + Quick_Cross + Average_Movement + Back_n_forth + Trace + unknown + Trapped), the total number of altered encounters (\code{alt_enc} = Bounce + Back_n_forth + Trace + Trapped) and/or total number of individuals affected (\code{unique_ind}).
  Result of this expression is transformed to range from 0 to 1 using equation: 
  \deqn{function(x){(x-min(x))/(max(x)-min(x))}}
  Default is [scenario in paper....]}
  \item{show_plot}{\code{logical}, \code{TRUE} if the barrier ranking should be visualised in the ploting device.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an object of class \code{sf} and \code{data.frame} containing one row per fence ID, and columns with information on:
\item{unique_ind}{number of individuals that encountered the barrier.}
\item{Quick_Cross, Bounce, unknown,
Avergae_movement, Trace, Back_n_forth,
Trapped}{total number of the different event types.}
\item{total_enc}{total number of encounters.}
\item{alt_enc}{total number of altering-movement encoutners (\code{Bounce + Back_n_forth + Trace + Trapped}).}
\item{index}{ranking index as calculated by \code{index_fun} and transformed to a scale of 0 to 1.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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

results_pron <- BaBA(animal = pronghorn, barrier = fences, interval = 2, b_hours = 4, p_hours = 36, d = 110, max_cross = 4, exclude_buffer = T, w = 7, img_suffix = "PRON")

fence_ranking <- BaRanking(classification = results_pron$classification, barrier = fences, d = 110, min_total_enc = 10, index_fun = expression((alt_enc/total_enc)*unique_ind))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
