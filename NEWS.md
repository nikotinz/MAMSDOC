## Version 3.0.3
- Fixed DTL cumulative sample size calculation in summary function
- Fixed print output for summary for DTL
- Fixed subscript out of bounds error in separate stopping rules for J>K
- Changed the way to define the probability of 'T1 being best' when 'method = "sep"', now matching the logic of other methods.

## Version 3.0.2
- Fixed Cohen's d effect size estimates in the summary() function for "MAMS" objects
- Fixed sample size calculation when J=1 and effect sizes were specified via 'delta' arguments (i.e., delta, delta0 and sd)

## Version 3.0.1
- Added handling allocation ratios at the first stage when the control arm is greater than 1 (thanks to Peter Greenstreet for the feedback)
- Added extra checks for input parameters of allocation ratios
- Fixed 'maximum iteration number' computation
- Renamed rows and cols in rMat for 'mams' object

## Version 3.0.0
- Added 2 new designs:
  - method = 'separate' for MAMS with a separate stopping rule
  - method = 'dtl' for drop-the-losers design
  - method = 'simultaneous' (default) corresponds to the MAMS with a simultaneous stopping rule
- Refined the summary(), print() and plot() functions for "MAMS" objects accordingly
- Updated the C code to compile with STRICT_R_HEADERS=1 (mandatory for R >= 4.5.0)

## Version 2.0.2
- NA related C code bug fix related to gcc compilers
- Add documentation for better clarity and user guidance

## Version 2.0.1
- Bug fix when K=1 in mams.sim
- Add src/init.c file

## Version 2.0.0
- Parallelisation framework included for several functions
- Recoding of parts of the functions in C
- Inclusion of new package help file
- Update of authors

## Version 1.4.2
- Fixed an error in tite.mams
- Update of author contact information

## Version 1.4.1
- Minor bug fix

## Version 1.4
- Adjusted code to also provide answers for K=1 and J=1
- Fixed an error in tite.mams
- Added warning to show that sample size search stopped due to maximum reached

## Version 1.3
- References updated

## Version 1.2
- New function 'ordinal.mams' for ordinal and binary endpoints
- New function 'tite.mams' for time-to-event endpoints
- Effect sizes can be specified on traditional or probability scale
- Search for sample size is capped at a maximum value
- Harmonised function and input names (most notably, 'step_down_mams' and 'update_mams' are now called 'stepdown.mams' and 'stepdown.update')
- Corrected minor bugs in 'mams'

## Version 1.1
- Corrected minor bug in mams.sim for K=1

## Version 1.0
- Altered mams function to improve consistency for different allocation ratios

## Version 0.9
- Corrected a bug in update_mams

## Version 0.8
- Removal of some depreciated code

## Version 0.7
- Modified print and summary to only return integer sample sizes
- Updated mams function to deal with minor inconsistency when J=1

## Version 0.6
- Corrected a bug in new.bounds
- Corrected a bug in mams.sim
- A few minor bugs in the plot functions corrected
- Arguments bty and las for the plot functions added
- Inconsistencies in the documentation removed

## Version 0.5
- A bug in function mams.sim corrected

## Version 0.4
- A bug in function mams corrected

## Version 0.3
- New function (step_down_mams) to calculate stopping boundaries for all intersection hypothesis tests in a closed testing procedure.
- New function (update_mams) to update boundaries at an interim analysis to take account of unplanned treatment selection and/or sample size reestimation.
- Documentation updated

## Version 0.2
- New function (mams.sim) to simulate mams studies included
- New function (new.bounds) to update boundaries based on observed number of observations included
- Citation information updated
- Output data.frame of function mams generalized
- Additional option to suppress sample size calculation included in function mams
- Documentation updated

## Version 0.1
- Initial release
