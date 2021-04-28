\name{BaBA}
\alias{BaBA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Barrier Behavior Analysis (BaBA)
}
\description{
This function classifies animal movement near linear barriers into 6 different behavioral categories: quick cross, average movement, bounce, back-and-forth, trace, and trapped. The classification unit of BaBA first identifies \emph{'encounter events'}, defined by a series of continuous locations that fall within the barrier buffers. Each event will then be classified into one of the 6 barrier behaviors. Even though many arguments have default settings, we recommend users adjust these values to best fit the ecology and management goals of the studied species. 
}
\usage{
BaBA(animal, barrier, d,
interval = NULL, b_time = 4, p_time = 36, w = 168, tolerance = 0, units = c("secs", "mins", "hours", "days", "weeks"),
max_cross = 0,  sd_multiplier = 1, exclude_buffer = F,
export_images = F, img_path = "event_imgs", img_suffix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{animal}{
A \code{"SpatialPointsDataFrame"} of animal movement with a column named \code{"Animal.ID"} that identifies each individual; a column named \code{"date"} in \code{"POSIXct"} format. \code{CRS} has to be the same as \code{barrier}, preferably in projected UTM coordinate system. More see Details.
}
\item{barrier}{
A \code{"SpatialLinesDataFrame"} showing barrier locations in the area overlapped with \code{animal} movement data. In the same spatial projection as \code{animal}.
}
\item{d}{
Barrier buffer size in meters if \code{barrier} has a projected coordinate system CRS, in the units of \code{barrier} map otherwise. 
}
\item{interval}{
Time interval of the movement data (unit specified in \code{units}). If not specified, BaBA will use the most frequent time difference between steps as \code{interval}, which might affect result accuracy if the data has irregular time interval or a lot of missing data. 
}
\item{b_time}{
Maximum duration, in the same unit as \code{interval}, that an encounter event would be considered as a short event (\emph{'bounce'} or \emph{'quick cross'}). Must be divisible by \code{interval}

}
\item{p_time}{
Minimum duration, in the same unit as \code{interval}, that an encounter event would be considered as a \emph{'trapped'} condition. Must be divisible by \code{interval}.
}
\item{w}{
The length of time, in the same unit as \code{interval}, to include around the encounter event to calculate average movement straightness using a moving window. Locations included are all locations within \code{w/2} before the first location of the encounter event and \code{w/2} after the last location of the event. More see Details.
}
\item{tolerance}{
The maximum duration, in the same unit as \code{interval}, to allow points that are outside of barrier buffer but between 2 sets of points within the buffer to be included, so that all points are considered as one continous \emph{encounter event}. Useful when movement data has a very high temporal resolution. Must be divisible by \code{interval}.
}
\item{units} {
The temporal units of \code{interval}, \code{b_time}, \code{p_time}, \code{w}, and \code{tolerance}. One of "secs", "mins", "hours", "days", and "weeks".
}
\item{max_cross}{
The maximum number of crosses in an encounter event allowed for in tracing and back-and-forth behavior. More see Details. 
}
\item{sd_multiplier}{
A \code{numeric} value to determine the numbers of standard deviations used that defines the normal range of movement straightness. Default is 1. More in Details. 
}
\item{exclude_buffer}{
\code{Logical}. Whether to consider movement locations within barrier buffers when calculating average movement straightness. More see Details.  
}
\item{export_images}{
\code{Logical}. If \code{TRUE,} will export snapshot of event locations as \code{.png} files named with the classification of the event, the \code{LocationID} of the animal, the \code{burstID} of the event and potential suffix specified by argument \code{export_suffix}.
}
\item{export_path}{
When \code{export_images} is \code{TRUE}, the \code{character} name of the folder that the images should be exported to. The folder will be created in the working directory if it does not exist.
}
\item{img_suffix}{
\code{character} to be added at the end of the naming if the images are exported when \code{export_images} is \code{TRUE}
}
}
\details{
\code{BaBA} works better with cleaned \code{animal} dataframe. Make sure your date column is cleaned so that the time difference between any step is divisible by \code{interval}. It is OK to have missing steps but remove repetitive timesteps before running \code{BaBA}. 

Barrier buffer distance \code{d} represents the distance on both sides of the barrier lines within which animal movement locations are considered 'encounters' with the barrier, and continuous locations form an \emph{'encounter event'}. This might be the most important parameter to set the BaBA. The buffer distance will affect numbers and durations of trajectories identified as 'encounter events'. For species with different movement capacities or ecological attributes, barrier effect distance might be different. You can decide the distance either by a priori knowledge or testing a range of distance and compare the results.

An average straightness is used to be compared to the straightness of the straightness of a to-be-classified encounter event in order to determine whether the event is \emph{'trace'} (more straightness than average), \emph{'back-n-forth'} (less straight than average), or \emph{'normal movement'} (similar to average). The reason to apply a moving window method with a width of \code{w} to calculate the average straightness is because some animal movements, such as migrations, show great seasonal variations. Moving window method calculates localized average straightness around the time when the to-be-classified encounter event occurs. 

The default average straightness calculation considers movement segments \code{w/2} days before and after the focal encoutner event, which can include movement locations within and outside of the barrier buffers. One might want to control potential impacts of barriers on movement and only calculate average movement straightness based on locatons outside of fence buffers. Use \code{exclude_buffer} to indicate whether or not to include movement locations inside barrier buffers. The default is \code{FALSE} because in the sample study fence density is relatively high and animals fall in barrier buffers all the time. If excluded, not enough continuous locations within the time window would be outside of the buffer.

The average and standard deviation of straightness measurements calculated by the moving window method are used to define the normal range of movement straightness. Any encounter event with a straightness < (average straightness + \code{sd_multiplier} * sd straightness) would be classified as \emph{'trace'} (more straight than normal), < (average straightness - \code{sd_multiplier} * sd straightness) would be \emph{'back-and-forth'} (more tortuous than normal), and in between would be \emph{'average movement'}. To make the standard more strict for an event to be classified as not 'normal' (i.e. \emph{'back-and-forth'} or \emph{'trace'}), use a larger number for the \code{sd_multiplier}.

When the barriers are curvy or the temporal resolution of the movement data is coarse, straight lines between movement locations (movement segments) might appear to be crossing the barrier even if the locations are on the same side of the barrier. \code{max_cross} allows some numbers of intersections between movement segments and barriers to be included in the \emph{'trace'} and \emph{'back-and-forth'} behavior. When the intersections are larger than \code{max_cross}, the encounter event will be classified as \emph{'unknown'}. 

}
\value{
\code{BaBA} returns a list.
\item{$classification}{
a dataframe. Each row represents one encounter event. \code{AnimalID} is the ID of the individual; \code{burstID} is the ID of the encounter event and it matches the ID of the exported images; \code{easting} and \code{northing} are the coordinates of the starting location of the encounter event; \code{start_time}, \code{end_time}, and \code{duration} describe the temporal characteristics of the encounter event; \code{cross} is the number of intersections between the encounter event trajectory and the barrier; \code{straightness} is the straightness of the encounter event trajectory and is NaN when the encounter event is \emph{'bounce} or \emph{'quick cross'}; finally, \code{eventTYPE} is the barrier behavior classification of the encounter event.
}
\item{$encounters}{
a \code{SpatialPointsDataFrame} showing the locations of classified encounter events represented by the starting location.
}
}
\references{
Xu W, Dejid N, Herrmann V, Sawyer H, Middleton AD. Barrier Behaviour Analysis (BaBA) reveals extensive effects of fencing on wide-ranging ungulates. J Appl Ecol. https://doi.org/10.1111/1365-2664.13806
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
library(BaBA)

# load sample data.
data("pronghorn")
data("fences")

# Animal data must be a "SpatialPointsDataFrame" created by the "sp" package. 
# individual ID column should be named "Animal.ID" and timestamp column should be named "date"
class(pronghorn)
names(pronghorn)

# run BaBA on the pronghorn data
results_prong <- BaBA(animal = pronghorn, barrier = fences, d = 110, max_cross = 4)

# View BaBA results
View(results_prong$classification)

# plot encounter event locations
plot(fences)
plot(results_prong$encounters, add = T)

# write the encounters as shapefile
writeOGR(results_prong$encounters, ".", "encountersPRONG", driver = "ESRI Shapefile", overwrite_layer = T)

# export event images to visually check event classifications. Using mule deer data as an example
data("muleDeer")
results_deer <- BaBA(animal = muleDeer, barrier = fences, d = 90, export_images = T, img_suffix = "DEER")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
