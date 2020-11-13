# MA-of-continuous-outcomes-with-baseline-imbalance

## Description 
In an parallel-group RCT, continuous outcome data are typically measured at two timepoints; pre- and post-treatment. When synthesizing trials of these data, the effect size of interest is commonly the mean difference, which can be calculated under various methods.

We provide code to fit standard aggregate data approaches: follow-up scores and change scores and less standard approaches (originally implemented in STATA), i.e., the [Trowman method](10.1016/j.jclinepi.2007.03.014) and [modified Trowman method](https://doi.org/10.1002/sim.5726). The less standard approaches are not recommended.

In addition, we propose the preferred approach of generating pseudo IPD from the appropriate summary data (*sufficient statistics*) and applying ANCOVA in an one-stage or two-stage model.


## Codes 
* [R](https://github.com/Katerina-Pap/MA-of-continuous-outcomes-with-baseline-imbalance/blob/main/R/pseudoIPD.pre-post.R)
* [SAS](https://github.com/Katerina-Pap/MA-of-continuous-outcomes-with-baseline-imbalance/blob/main/SAS/pseudoIPD.pre-post.sas)
* [STATA](https://github.com/Katerina-Pap/MA-of-continuous-outcomes-with-baseline-imbalance/blob/main/STATA/pseudoIPD.pre-post.do)
