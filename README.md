# Supporting information for "Allowing growing dimensional binary outcomes via the multivariate probit Indian buffet process"

Code and data to implement the analysis in [Allowing growing dimensional binary outcomes via the multivariate probit Indian buffet process.](https://arxiv.org/abs/2402.13384)
The repository contains:
1. The file `MVPIBPh_nocov.cpp`implements the two-stage algorithm to estimate the hierarchical MVP-IBP model, and `MVPIBPh.cpp` implements the corresponding extension to covariate-dependent model. The file `gibbs_factorMVPIBP.R` implements the Gibbs sampler to estimate the factor MVP-IBP.
2. The file `IntroPlots.R` contains the code to reproduce Figure 1 and Figure 2 in the manuscript. 
3. To run the simulation studies, see `SimStudy.R`.
4. The application to the fungi data is implemented in `fungiAnalysis.R`, while to replicate the prediction results for p* shown in Table 1 see `prediction_pstarn.R`. 
