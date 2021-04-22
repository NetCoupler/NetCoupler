
# Notes from meeting

- Instead of making sim data, just Fabian can install on server

- compare implementation of below (Clemens version) vs the current (my) version
    - From Clemens: I have checked in the ceramides project– the criteria for direct effect identification were i) P < 0.05 in the non-neighbor-adjusted model (rule1), ii) P < 0.1 in all neighbor-adjusted models (rule 2), iii) same sign of the beta in all the models (rule 3). I did not correct for multiple testing in the primary selection process.
    It would be interesting to apply the current NetCoupler-version to the Ceramides project with these same selection cutoffs and reproduce the results. In principle, it would be good to make these cutoffs' manual setting amendable as an option.

- As idea: provide argument for options as list.

- Determine journal and decide on whether it is method/theory vs application paper

- Schedule next meeting
    - agenda etc
    
- Surv can have three arguments
    - Or do it as a pre-processing (Surv or strata in mutate)

- things throw off when adjusting for variables with similar names (e.g "cer" when there is a "cer14.0" in the network)
    - based on some pattern selection?
    
- coxph cluster() within the formula interface? Same with strata (but that could be put in preprocessing)
    - cluster doesn't work.
