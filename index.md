<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/MAMS)](https://cran.r-project.org/web/packages/MAMS)
[![Downloads](https://cranlogs.r-pkg.org/badges/MAMS?color=blue)](https://cran.rstudio.com/package=MAMS)
<!-- badges: end -->

# MAMS <img src="./vignettes/logo.png" width="100" align="right" alt="MAMS logo"/>

<h4>Designing multi-arm multi stage studies</h4>

</span>

<br>

This package allows to design multi-arm multi-stage (MAMS) studies with
asymptotically normal endpoints and known variance. It considers normal, binary,
ordinal and time-to-event endpoints in which either the single best treatment or
all promising treatments are continued at the interim analyses.

### Installation

You can install the latest released version from
[CRAN](https://cran.r-project.org/package=MAMS) from within R:

```r
install.packages("MAMS")
```

### Details

Currently implemented functions are:

-   **`mams()`**: a function allowing to design multi-arm multi-stage studies with
    normal endpoints,

-   **`new.bounds()`**: a function allowing to update the lower and upper boundaries
    of a multi-arm multi-stage study, typically initally defined by `mams()`,
    based on observed sample sizes,

-   **`mams.sim()`**: a function allowing to simulate multi-arm multi-stage studies
    given chosen boundaries and sample size, and estimates power and expected
    sample size,

-   **`stepdown.mams()`**: a function allowing to find stopping boundaries for a 2- or
    3-stage (stepdown) multiple-comparisons-with-control test,

-   **`stepdown.update()`**: a function allowing to update the stopping boundaries of
    a multi-arm multi-stage study, typically initally defined by
    `stepdown.mams()`, at an interim analysis as well as allowing for unplanned
    treatment selection and/or sample-size reassessment,

-   **`ordinal.mams()`**: a function allowing to design multi-arm multi-stage studies
    with ordinal or binary endpoints,

-   **`tite.mams()`**: a function allowing to design multi-arm multi-stage studies
    with time-to-event endpoints.  
    
    We refer to Jaki et al (2019) for an overview
    of the package as well as to Magirr et al (2012) and Magirr et al (2014) for
    theoretical details.


