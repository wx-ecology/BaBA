# BaBA
Public Script for Barrier Behavior Analysis (BaBA). This analysis is to detect and quantify animal movement behavior upon encountering with linear barriers (e.g. roads, fences, pipelines, rivers).

The 7 kinds of behaviors that BaBA classifies is: Bounce, Trace, back-n-forth, Average Movement, Quick Cross, Trapped, and unknown. See illustration:
![BaBA catogory classes](https://github.com/wx-ecology/BaBA/blob/master/BaBA_Catogories.png)

To conduct this analysis, you need animal tracking data, then follow the next two script. 

*01-Generate-straightness-table*: The differentiations among Trace, average movement, and back-n-forth are based on movement trajectory straightness. This script is to generate a straightness table that calculate trajectorys straightness with duration from b (default set to 6) hours to p hours (default to 36 hours) using a moving window method. With this table, avergae movement straghtness of certain time can be calculated and straightness of barrier encountering event will be compared with this weekly average straightness. 

*02-Barrier-Behavior-Classification*: this script is to identify encounter events and to classify them into one of the seven classes. The output is a table with animal ID, encounter event ID (time of event), event coordinates, duration, # of crosses, and event type.

Relevant publication: 
In prep. 

Last update: Dec 19, 2019 
Wenjing Xu
