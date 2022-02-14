# InSilicoVA
[![R-CMD-check](https://github.com/verbal-autopsy-software/InSilicoVA/workflows/R-CMD-check/badge.svg)](https://github.com/verbal-autopsy-software/InSilicoVA/actions) [![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/InSilicoVA)](https://cran.r-project.org/package=InSilicoVA) [![](https://cranlogs.r-pkg.org/badges/InSilicoVA)](https://cran.r-project.org/package=InSilicoVA) [![](https://cranlogs.r-pkg.org/badges/grand-total/InSilicoVA?color=orange)](https://cran.r-project.org/package=InSilicoVA)

To install the package from CRAN, run the script: 
```r
install.packages("InSilicoVA")
```

To use InSilicoVA, 'rJava' package is required and it requires Java 6 or above properly installed and linked to R. See <a href="Documents/Insilico-manual.pdf">the user manual</a> for more detail.

### Reference
Tyler H. McCormick, Zehang Richard Li, Clara Calvert, Amelia C. Crampin, Kathleen Kahn, and Samuel J. Clark. <a href="http://arxiv.org/abs/1411.3042">Probabilistic Cause-of-death Assignment using Verbal Autopsies</a>. _Journal of the American Statistical Association_ 111.515 (2016): 1036-1049.



### What's New
- [CRAN version](https://cran.r-project.org/web/packages/InSilicoVA/news/news.html) 
- [Current developer version](InSilicoVA/NEWS.md)
 

### Additional information when running large models
If you hit the error message compaining about java heap size, e.g. ``java.lang.OutOfMemoryError: Java heap space``

Try start a new R session and run this line of code **before** calling ``library(openVA)`` or ``library(InSilicoVA)``:

> options(java.parameters = c("-XX:+UseConcMarkSweepGC"))
