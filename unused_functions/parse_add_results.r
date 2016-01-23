Parse_Add_Result <- function(N_sub.j, C.j, S.j, N_level.j, pool.j, fit, fit.add){
###########################################################
# function to add Java output to previous results and parse into correct place 
# @param:
#	   various java arguments
#	   both new and old Java output
# @values:
#		list of variables parsed 
    # remove the first row from the added matrix, since it is the starting point
    fit.all <- rbind(fit, fit.add[-1, ])
    # parse all the saved iterations
    out_all <- ParseResult(N_sub.j, C.j, S.j, N_level.j, pool.j, fit.all)
    # parse only the added fit to find the most recent values
    out_last <- ParseResult(N_sub.j, C.j, S.j, N_level.j, pool.j, fit.add)
    # replace the most recent values with the correct set
    out_all$mu.last <- out_last$mu.last
    out_all$sigma2.last <- out_last$sigma2.last
    out_all$theta.last <- out_last$theta.last     
    return(out_all)
}