# bootstep

A Stata program to perform bootstrap estimation for two-step estimations with a generated regressor used as independent variable"

by David Veenman (University of Amsterdam)

`bootstep` provides a procedure to obtain bootstrapped standard errors when the predicted values from a first-step estimation are used as independent variable in a second-step estimation. A shown in prior work, ignoring the sampling error of the generated regressor causes the second-step standard errors to be understated in the second-step estimation (see [Pagan 1984](https://doi.org/10.2307/2648877) and [Murphy and Topel 1985](https://doi.org/10.1198/073500102753410417); see [Chen et al. (2020)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3724730) for illustrations in accounting/finance research). The procedure allows for estimation of standard errors that are robust to clustering in one or two dimensions using a pairs-cluster bootstrap approach and the formula for twoway cluster-robust standard errors from [Thompson (2011)](https://doi.org/10.1016/j.jfineco.2010.08.016) and [Cameron, Gelbach, and Miller (2011)](https://doi.org/10.1198/jbes.2010.07136). 

---

The program can be installed by simply saving the \*.ado file into your local Stata directory that contains additional ado programs. To identify this folder, type and execute "sysdir" in your Stata command window and go to the folder listed at "PLUS:". Make sure you place the program in the folder with the correct starting letter of the program (i.e., the folder named "b") and the file extension is correctly specified as \*.ado.

---

Version history and changes:

  2021dec03: version 1.0.0

---

Program options:

- *nboot* (required): number of bootstrap replications.  
- *cluster(varlist)* (optional): allows for bootstrap estimation of cluster-robust standard errors in one or two dimensions. 
- *logit* (optional): allows the first step to be estimated using logit (provided the first-step dependent variable is an indicator variable); the first-step fitted values are transformed into predicted probabilities.
- *absorb* (optional): allows for one dimension of fixed effects to be absorbed in the second-step estimation (in this case, the program relies on Stata's `areg`).
- *standardize* (optional): standardizes the coefficient on the generated regressor in the second-step estimation to allow for a direct comparison with the case in which the first-step dependent variable would be included directly as an independent variable in the second-step estimation. The standardization follows from equation (3) of [Hou and Loh (2016)](https://doi.org/10.1016/j.jfineco.2016.02.013).

---

Examples of output:


