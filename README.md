# InSilicoVA  [![Build Status](https://travis-ci.org/richardli/InSilicoVA.svg?branch=master)](https://travis-ci.org/richardli/InSilicoVA)

To install the package from CRAN, run the script: 
```r
install.packages("InSilicoVA")
```

To use InSilicoVA, 'rJava' package is required and it requires Java 6 or above properly installed and linked to R. See <a href="Documents/Insilico-manual.pdf">the user manual</a> for more detail.

### Reference
Tyler McCormick, Zehang Li, Clara Calvert, Amelia Crampin, Kathleen Kahn, and Samuel Clark. <a href="http://arxiv.org/abs/1411.3042">Probabilistic Cause-of-death Assignment using Verbal Autopsies</a>. <em>To appear, Journal of the American Statistical Association</em>. 



## What's New
#### version 1.1.5 (2017-11-17)
- Fix issue for data without external causes.
- Fix issue with print methods after last update.
- Fix issue with neonate and child death assigned to impossible external causes.
- Change data check steps to update symptoms to be missing instead of absence to introduce symptom dependence structures. 
- Automatic remove impossible causes from CSMF based on subpopulation gender and age. 

#### version 1.1.4 (2017-01-24)
- Fix issue with sub-population specification and output labels of non-standard InterVA-4 input.

#### version 1.1.3 (2017-01-02)
- Fix issue with label order in probbase output.
- Allow input "Yes" to be denoted with either "y" or "Y".
- Default data checking rules to be consistent with InterVA-4.03 instead of the previous version.
- Fix typo in document.

#### version 1.1.2
- Bug fix for function csmf.diag()

#### version 1.1.1
- New function: updateIndiv(): Now the C.I. for individual COD distribution is not calculated by default, in order to save computation time. But instead it is updated using this new function. 
- Minor updates of function syntax 
- Improved reliability for larger datasets 
- Update version numbering convention
