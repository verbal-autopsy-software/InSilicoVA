# InSilicoVA 

To install the package from CRAN, run the script: 
```r
install.packages("InSilicoVA")
```

To use InSilicoVA, 'rJava' package is required and it requires Java 6 or above properly installed and linked to R. See <a href="Documents/Insilico-manual.pdf">the user manual</a> for more detail.

### Reference
Tyler McCormick, Zehang Li, Clara Calvert, Amelia Crampin, Kathleen Kahn, and Samuel Clark. <a href="http://arxiv.org/abs/1411.3042">Probabilistic Cause-of-death Assignment using Verbal Autopsies</a>. <em>To appear, Journal of the American Statistical Association</em>. 



## What's New

#### version 1.1.1
- New function: updateIndiv(): Now the C.I. for individual COD distribution is not calculated by default, in order to save computation time. But instead it is updated using this new function. 
- Minor updates of function syntax 
- Improved reliability for larger datasets 
- Update version numbering convention
