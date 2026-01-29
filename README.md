The purpose of this code is to replicate the analyses associated with the manuscript:

Parvin, L.K., V. Jirinec, M. Leu, and J.J. Valente. (in review). Post hoc statistical analyses cannot standardize occupancy estimates for mobile species across variable survey protocols. _Ecological Solutions and Evidence_.__

This manuscript simulates point count sampling of a simulated Wood Thrush (_Hylocichla mustelina_)__ population and evaluates whether asymptotic regression models can be used to recover true instantaneous and daily occupancy rates. The script for simulating Wood Thrush locations and the resulting point count data is already published and available via Appendix 5 at (https://doi.org/10.5066/P98BJAQU). That script produces csv files named simResults# where the # sign refers to the simulation iteration. For ease of use, we have uploaded results from two simulation iterations here as CSV files named "simResults1.csv" and "simResults2.csv".

The first step in replicating the analyses associated with Parvin et al. (in review) involves opening "pointCountOccupancyAnalysis.R"

asymptoticModelsTmp

asymptoticModeling.R must be run three times, once for each density (can be specified at the top of the script). binomialRegression.R should be run subsequently.
===

