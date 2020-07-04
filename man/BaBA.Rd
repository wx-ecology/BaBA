\name{BaBA}
\alias{BaBA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Barrier Behavior Analysis (BaBA)
}
\description{
This function classifies animal movement near linear barriers into 6 different behavioral categories: quick cross, average movement, bounce, back-and-forth, trace, and trapped.
}
\usage{
BaBA(animal, barrier, d = 50,
interval = 2, b_hours = 4, p_hours = 36,
w = 7, max_cross = 4, tolerance = 0, sd_multiplier = 1, exclude_buffer = T,
 export_images = T,img_path = "event_imgs", img_suffix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{animal}{
A \code{SpatialPointsDataFrame} of animal movement with a field that identifies each indvidual: \code{Animal.ID}, a \code{date} column in \code{POSIXct} format. \code{CRS} has to be the same as \code{barrier}, preferably in projected coordinate system.
}
\item{barrier}{
\code{SpatialLinesDataFrame} showing barrier locations in the area overlapped with \code{animal} movement data. In the same spatial projection as \code{animal}.
}
\item{d}{Barrier buffer size in meters if  \code{barrier} has a longitude/latitude CRS, in the units of \code{barrier} map otherwise. 
}
\item{interval}{
Time interval (in hours) of the movement data.
}
\item{b_hours}{
Maximum encounter duration, in hours, that an event would be considered as a 'short encounter' ('Bounce' or 'Quick Cross'). Must be divisible by \code{interval}

}
\item{p_hours}{
Minimum encounter duration, in hours, that an event would be considered as a "trapped" condition. Must be divisible by \code{interval}.
}
\item{w}{
The number of days to include around the encounter event to calculate average movement straightness using a moving window. Locations included are all locations within \code{w/2} days before the first location of the encounter event and \code{w/2} days after the last location of the event.}
\item{max_cross}{
The maximum number of crosses allowed for in tracing and back-and-forth behavior. Sometimes, depending on the shape of the barrier and temporal interval of movement data, the straight segments between movement locations might appear to be crossing the barrier even if the locations are actually on the same side of the barrier. This parameter allows some crossings to be included in the 'trace' behavior.
}
\item{tolerance}{
The maximum duration, in hours, to allow points that are outside of barrier buffer but between 2 sets of points within the buffer to be included as one long event. Must be divisible by \code{interval}
}
\item{sd_multiplier}{
A \code{numeric} value to determine the range of standard deviations to consider a movement straightness within, below, or beyond the normal range. The average and standard deviation of straightness measurements is calculated by the moving window method. Any encounter event with a straightness < than average straightness + \code{sd_multiplier} * sd straightness will be classified as Trace (more straight than normal), < than average straightness - \code{sd_multiplier} * sd straightness it would be Back-and-forth (more tortuous than normal), and in between would be Average movement. 
}
\item{exclude_buffer}{
\code{Logical}. Whether to consider movement locations within barrier buffers when calculating average movement straightness.  
}
\item{export_images}{
\code{logical}. If \code{TRUE,} will export snapshot of event locations as \code{.png} files named with the classification of the event, the \code{LocationID} of the animal, the \code{burstID} of the event and and potential suffix specified by argument \code{export_suffix}.
}
\item{export_path}{
When \code{export_images} is \code{TRUE}, the \code{character} name of the folder the images should be exported to. The folder will be created in the working directory if it does not exist.
}
\item{img_suffix}{
\code{character} to be added at the end of the naming if the images exorted when \code{export_images} is \code{TRUE}
}
}
\details{
Barrier buffer distance \code{d} represents the distance on both sides of the barrier lines within which animal movement locations are considered "encounters" with the barrier. This might be the most important parameter to set the BaBA. The buffer distance will affect numbers and durations of trajectories identified as 'encounter events'. For species with different movement capacities or ecological attributes, barrier effect distance might be different. You can decide the distance either by a priori knowledge or testing a range of distance and compare the results.

An average straightness is used to be compared to the straightness of the straightness of a to-be-classified encounter event in order to determine whether the event is 'trace' (more straightness than average), 'back-n-forth' (less straight than average), or 'normal movement' (similar to average). The reason to apply a moving window method with a width of \code{w} to calculate the average straightness is because some animal movements, such as migrations, show great seasonal variations. Moving window method calculates localized average straightness around the time when the to-be-classified encounter event occurs. The default average straightness calculation considers movement segments \code{w/2} days before and after the focal encoutner event, which can include movement locations within and outside of barrier buffers. One might consider movements within the buffer to have been affected by the barrier and will bias what are considered as "average movement". Use \code{exclude_buffer} to indicate whether or not to include movement locations inside barrier buffers. The default is \code{FALSE} because in the sample data fence density is relatively high and animals fall in barrier buffers all the time. If excluded, not enough continuous locations within the time window would be outside of the buffer.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Xu, W., Nandintsetsegm D., Sawyer, H., Herrmann, V., Middleton, A. (n.d.). Barrier Behavior Analysis (BaBA) reveals extensive effects of fencing on wide-ranging animals
}
\author{
%%  ~~authors~~
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
data("muleDeer")
data("fences")


results_prong <- BaBA(animal = pronghorn, barrier = fences, interval = 2, b_hours = 4, p_hours = 36, d = 110, max_cross = 4, exclude_buffer = T, w = 7, img_suffix = "PRONG")
results_deer <- BaBA(animal = muleDeer, barrier = fences, interval = 2, b_hours = 4, p_hours = 36, d = 110, max_cross = 4, exclude_buffer = T, w = 7, img_suffix = "DEER")

View(results_prong$classification)
View(results_prong$encounters@data)
plot(results_prong$encounters)

View(results_deer$classification)
View(results_deer$encounters@data)
plot(results_deer$encounters)

# write the encounters as shapefile
writeOGR(results_prong$encounters, ".", "encountersPRONG", driver = "ESRI Shapefile", overwrite_layer = T)
writeOGR(results_deer$encounters, ".", "encountersDEER", driver = "ESRI Shapefile", overwrite_layer = T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
