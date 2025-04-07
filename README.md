# Barrier Behaviour Analysis

![Version](https://img.shields.io/badge/version-2.3-blue)

## Description
Barrier Behaviour Analysis (BaBA) is a spatial- and temporal-explicit method to identify and classify barrier behaviors, fence behaviors in our case, based on GPS tracking data and linear spatial features. Barrier behaviors can be used to examine permeability of barriers for animal movement.

The 6 kinds of behaviors that BaBA classifies is: bounce, trace, back-n-forth, average movement, quick cross, trapped. See illustration:
![BaBA catogory classes](BaBA_Catogories.png)

To install the latest development version of BaBA, in an R session, type: 

**devtools::install_github("wx-ecology/BaBA")**

## BaBA Workflow
![BaBA workflow](Flowchart.png)

## Updates
V2.3 (Apr 7, 2025)
1. fixed 'fetched pts' issue and ID mismatch (thank you A. Seveque!). 

V2.2 (Dec 5, 2024)
1. improved visuals of exported images for "bounce" and "quick cross".

V2.1 (Apr 23, 2024)
1. fixed issues related to straightness calculation that produces excessive "unknown" results
2. fixed issues related to tolerance in the new sf-based function 

V2.0 (Sep 27, 2023) 
1. Update package to run on sf instead of sp
2. details [here](https://github.com/wx-ecology/BaBA/pull/4#issue-1903676738)

V1.3 (Apr 14, 2023)
1. Fixed bug related to "tolerance"

V1.2 (May 26, 2022)
1. Solve encounters and classified encounters mismatch issue

V1.1 (Apr 10, 2021)
1. Improved the visualization of exported images 
2. Added time unit option and unified all temporal parameter units

## Relevant publication and example applications
Xu W, Dejid N, Herrmann V, Sawyer H, Middleton AD. (2021). Barrier Behaviour Analysis (BaBA) reveals extensive effects of fencing on wide-ranging ungulates. Journal of Applied Ecology. 58(4), 690-698. https://doi.org/10.1111/1365-2664.13806

Xu, W., Gigliotti, L. C., Royauté, R., Sawyer, H., & Middleton, A. D. (2023). Fencing amplifies individual differences in movement with implications on survival for two migratory ungulates. Journal of Animal Ecology, 92(3), 677-689. https://doi.org/10.1111/1365-2656.13879 

Aikens, E. O., Merkle, J. A., Xu, W., & Sawyer, H. (2025). Pronghorn movements and mortality during extreme weather highlight the critical importance of connectivity. Current Biology.https://doi.org/10.1016/j.cub.2025.03.010 
