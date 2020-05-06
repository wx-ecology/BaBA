\name{BaBA}
\alias{BaBA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Barrier Behavior Analysis (BaBA)
}
\description{
This function classifies animal movement near linear barriers into 7 different behavioral categories: Bounce, Quick Cross, Back-and-forth, Trace, Trap, Average Movement and  unknown.
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
\code{SpatialPointsDataFrame} of animal movement with a field that identifies each idvidual: \code{Location.ID}, a \code{date} column in \code{POSIXct} format. \code{CRS} has to be the same as \code{barrier}.
}
\item{barrier}{
\code{SpatialLinesDataFrame} showing barrier locations in the area overlapped with \code{animal} movement data. In the same spatial projection as \code{animal}.
}
\item{d}{Barrier buffer size in meters if  \code{barrier} has a longitude/latitude CRS, in the units of  \code{barrier} map otherwise. Distance on either side of the barrier lines within which animal movement is potentially affected by the presence of the barrier. This might be the most important parameter to set the BaBA. The buffer distance will affect numbers and durations of trajectories identified as 'encounter events'. For species with different movement capacities or ecological attributes, barrier effect distance might be different. You can decide the distance either by a priori knowledge, or testing a range of distance and compare the results.
}
\item{interval}{
Time interval (in hours) of the movement data.
}
\item{b_hours}{
Maximum encounter duration, in hours, that an event would be considered as a 'short encounter' ('Bounce' or 'Quick Cross'). Must be divisible by \code{interval}

}
\item{p_hours}{
Minimum encounter duration, in hours, that an event would be consiered as a "trapped" condition. Must be divisible by \code{interval}.
}
\item{w}{
Number of days to include around the encounter event  to calculate average movement straightness using a moving window. Locations included are all locations within w/2 days before the first location of the encounter event and w/2 days after the last location of the event. The average straightness is used to compare to the straightness of the focal encounter event in order to deterimine whether the event is 'trace' (more straightness than average), 'back-n-forth' (less straight than average), or 'normal movement' (similar to average). Note that the average straightness is calculated on all movement segments of the same duration of the encounter event within the time window set here, which include movement segments both within and outside of the barrier buffers. Use \code{exclude_buffer} to indicate whether or not to include movement locations inside barrier buffers.}
\item{max_cross}{Maximum number of crosses allowed for in tracing and back-and-forth behavior. Sometimes, depending on the shape of the barrier and temporal interval of movement data, the straight segments between movement locations might appear to be crossing the barrier even if the locations are actually on the same side of the barrier. This parameter allows some crossings to be included in the 'trace' behavior.
}
\item{tolerance}{
Maximum duration, in hours, to allow points that are oustide of barrier buffer but between 2 sets of points within the buffer to be included as one long event. Must be divisible by \code{interval}
}
\item{sd_multiplier}{
\code{numeric} value  used to multiply the standard deviation of all straightness measurements caluclated with the moving window method. Any encounter event with straightness > than average straightness + \code{sd_multiplier} * sd straightness will be classified as back-and-forth, < than average straightness - \code{sd_multiplier} * sd straightnessit would be Trace and inbetween would average movement. 
}
\item{exclude_buffer}{
\code{Logical}. When calculating average movement straightness, we consider movement segments \code{w} days before and \code{w} days after the focal encoutner event, which can include movement locations within and outside of barrier buffers. Should movement location within the buffer be excluded? default is \code{TRUE}, but you might consider changing this to \code{FALSE} if the barriers are relatively dense and animals fall in barrier buffers all the time and not enough continuous locations within the time window are outside of buffer (You would know that if you had a lot of 'unkown' classifications for events of medium durations.
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
Xu, W., Nandintsetsegm D., Sawyer, H., Herrmann, V., Middleton, A. (n.d.). Through the Devil's Rope: 
Barrier behavior analysis guides fencing modification prioritizations for wide-ranging ungulates
}
\author{
Wenjing Xu and Valentine Herrmann
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


results_pron <- BaBA(animal = pronghorn, barrier = fences, interval = 2, b_hours = 4, p_hours = 36, d = 110, max_cross = 4, exclude_buffer = T, w = 7, img_suffix = "PRON")


View(results_pron$classification)
View(results_pron$encounters@data)
plot(results_pron$encounters)

# write the encounters as shapefile
writeOGR(results_pron$encounters, ".", "encountersPRON", driver = "ESRI Shapefile", overwrite_layer = T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
