# Barrier Behavior Analysis
## Description
Barrier Behavior Analysis (BaBA)is a spatial- and temporal-explicit method to identify and classify barrier behaviors, fence behaviors in our case, based on GPS tracking data and linear spatial features. Barrier behaviors can be used to examine permeability of barriers for animal movement.

The 6 kinds of behaviors that BaBA classifies is: bounce, trace, back-n-forth, average movement, quick cross, trapped. See illustration:
![BaBA catogory classes](BaBA_Catogories.png)

To install the latest development version of BaBA, in an R session, type: 
**devtools::install_github("wx-ecology/BaBA")**

## BaBA Workflow
![BaBA workflow](Flowchart.png)

## Updates
V2.0 (September 27, 2023) 
1. Update package to run on sf instead of sp
2. details [here](https://github.com/wx-ecology/BaBA/pull/4#issue-1903676738)

V1.3 (Apr 14, 2023)
1. Fixed bug related to "tolerance"

V1.2 (May 26, 2022)
1. Solve encounters and classified encounters mismatch issue

V1.1 (April 10, 2021)
1. Improved the visualization of exported images 
2. Added time unit option and unified all temporal parameter units

## Relevant publication: 
Xu W, Dejid N, Herrmann V, Sawyer H, Middleton AD. Barrier Behaviour Analysis (BaBA) reveals extensive effects of fencing on wide-ranging ungulates. J Appl Ecol. 2020;00:1–9. https://doi.org/10.1111/1365-2664.13806

