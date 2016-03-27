# InSilicoVA-beta

Beta release of InSilicoVA package.

To install on mac, linux, and some Windows machine, run the script in R

```r
install.packages("ggplot2")
install.packages("coda")
install.packages("rJava")
install.packages("InSilicoVA_1.0.tar.gz", type = "source", repos = NULL)
```

For some Windows machine the source file might not work. Please let me know the detail (error message, your code, etc.) if that happens.

# What's New
#### version 1.1.1
- Minor updates of function syntax
- Improved reliability for larger datasets
- New function: updateIndiv(): Now the C.I. for individual COD distribution is not calculated by default, in order to save computation time. But instead it is updated using this new function.

