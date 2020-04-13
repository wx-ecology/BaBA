\name{BaBA}
\alias{BaBA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Barrier Behavior Analysis (BaBA)
}
\description{
This function classifies movement behaviors near barriers into 7 different categories: Bounce, Quick Cross, Back-and-forth, Trace, Trap, Average Movement and  unknown.
}
\usage{
BaBA(animal, barrier, d = 50,
interval = 2, b_hours = 4, p_hours = 36,
w = 7, max_cross = 4, tolerance = 0, exclude_buffer = T,
 export_images = T,img_path = "event_imgs", img_suffix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{animal}{
\code{SpatialPointsDataFrame} of animal movement with a field that identifies each idvidual: \code{Location.ID}, a \code{date} column in \code{POSIXct} format. \code{CRS} has to be the same as \code{barrier}.
}
\item{barrier}{
\code{SpatialLinesDataFrame} with the locations of all known barriers in the area covered by \code{animal} movement data. In the same spatial projection as \code{animal}.
}
\item{d}{ Barrier buffer in meters. Distance on either side of barrier lines within which an animal location/trajectory is considered to be potentially affected by the presence of the barrier. This might be the most important parameter to set the BaBA. The buffer distance will affect numbers and durations of trajectories identified as 'events'. For different species, barrier effect distance might be different. You can decide the distance either by experience, or testing a range of distance and compare the results.
}
\item{interval}{
Time interval (in hours) of the movement data.
}
\item{b_hours}{
Maximum encounter duration, in hours, that would be considered as a 'short encounter' ('Bounce' or 'Quick Cross'). Aka, what is quick for you? Must be divisible by \code{interval}

}
\item{p_hours}{
Minimum encounter duration, in hours, that you'd call the encounter event a "trapped" condition. Aka, what is too long for you? Must be divisible by \code{interval}.
}
\item{w}{
Number of days before the first location of the event and  and after the last one to consider when differenciating 'trace', 'back-n-forth' from 'normal movement'. We compare the straightness of the encounter event to the average straightness around the time the encounter event happens. How long would you set this window? Note that the average straightness is calculated on all movement segments of the same duration of the encounter event, within the window frame set here, befaore and after the encounter.}
\item{max_cross}{Maximum number of crosses allowed for in tracing behavior. Sometimes, depending on the shape of the barrier and the time interval, the straight segments between locations might appear to be crossing the barrier while the locations are actually always on the same side of the barrier. This parameter allows some crossing to happen in the 'trace' behavior.
}
\item{tolerance}{
not used.
}
\item{exclude_buffer}{
\code{Logical}. When considreing the average movement straightness looking at movement segments \code{w} days before and \code{w} days after the encoutner event, should locations within the buffer but from different encounter events be excluded? default is \code{TRUE} but you might consider changing this to \code{FALSE} if you have a very dense barriers and not enough continuous locations within the time window and outside of buffer (You would know that if you had a lot of 'unkown' classifications for events of medium durations.
}
\item{export_images}{
\code{logical}. If \code{TRUE,} will export snapshot of event locations as \code{.png} files named with the classification of the event, the \code{LocationID} of the animal, the \code{burstID} of the event and and potential suffix specified by argument \code{export_suffix}.
}
\item{export_path}{
When \code{export_images} is \code{TRUE}, the \code{character} name of the folder the images should be exported to. The folder will be created in the working directory if it does not exist.
}
\item{export_suffix}{
\code{character} to be added at the end of the naming if the images exorted when \code{export_images} is \code{TRUE}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
library(BaBA)
data("pronghorn")
data("fences")


results_pron <- BaBA(animal = pronghorn, barrier = fences, interval = 2, b_hours = 4, p.hours = 36, d = 110, max.cross = 4, exclude_buffer = T, w = 7, img_suffix = "PRON")


View(results_pron$classification)
View(results_pron$encounters@data)
plot(results_mule$encounters)

# write the encounters as shapefile
writeOGR(results_pron$encounters, ".", "encountersPRON", driver = "ESRI Shapefile", overwrite_layer = T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
