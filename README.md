# bootstep

A Stata program to perform bootstrap estimation for two-step estimations with a generated regressor used as independent variable

by David Veenman (University of Amsterdam)

**Preliminary: comments/feedback welcome**

`bootstep` provides a procedure to obtain bootstrapped standard errors when the predicted values from a first-step estimation are used as independent variable in a second-step estimation. A shown in prior work, ignoring the sampling error of the generated regressor causes the second-step standard errors to be understated in the second-step estimation (see [Pagan 1984](https://doi.org/10.2307/2648877) and [Murphy and Topel 1985](https://doi.org/10.1198/073500102753410417); see [Chen et al. (2020)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3724730) for illustrations in accounting/finance research). The procedure allows for estimation of standard errors that are robust to clustering in one or two dimensions using a pairs-cluster bootstrap approach and the formula for twoway cluster-robust standard errors from [Thompson (2011)](https://doi.org/10.1016/j.jfineco.2010.08.016) and [Cameron, Gelbach, and Miller (2011)](https://doi.org/10.1198/jbes.2010.07136). The program comes with a **fast** option based on mata optimization and [MacKinnon (2021)](http://qed.econ.queensu.ca/pub/faculty/mackinnon/working-papers/qed_wp_1465.pdf) that allows for much faster execution of the bootstrap resampling procedure. To limit the effects of downward bias in standard errors with a small number of clusters (G), the program applies a small sample correction to the variance matrix and relies critical values based on t(G-1) following the guidance of [Cameron and Miller (2015)](http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf).

---

The program can be installed by simply saving the \*.ado file into your local Stata directory that contains additional ado programs. To identify this folder, type and execute "sysdir" in your Stata command window and go to the folder listed at "PLUS:". Make sure you place the program in the folder with the correct starting letter of the program (i.e., the folder named "b") and the file extension is correctly specified as \*.ado.

---

Version history and changes:

  2021dec28 (version 1.0.1):
  - added option "fast" for much faster bootstrap execution using mata: the first-step regression estimates use the efficient computational algorithm of [MacKinnon (2021)](http://qed.econ.queensu.ca/pub/faculty/mackinnon/working-papers/qed_wp_1465.pdf), which speeds up the first-step estimations extremely; the second-step regression estimates cannot rely on this algorithm, but are also performed more efficiently given the use of mata functions.
  - improved inference by adding the small sample correction to the variance matrix and using critical values based on t(G-1) following [Cameron and Miller (2015)](http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf).

  2021dec03 (version 1.0.0): 
  - first commit to Github

---

Syntax:

**depvar1 indepvars1 | depvar2 indepvars2, options**

 - **depvar1** is the first-step dependent variable;
 - **indepvars1** is/are the first-step independent variable(s);
 - The vertical bar (**|**) separates the first and second step estimations;
 - **depvar2** is the second-step dependent variable;
 - **indepvars2** is/are the first-step independent variable(s).

---

Program options:

- **nboot** (required): number of bootstrap replications.  
- **cluster(varlist)** (optional): allows for bootstrap estimation of cluster-robust standard errors in one or two dimensions. 
- **logit** (optional): allows the first step to be estimated using logit (provided the first-step dependent variable is an indicator variable); the first-step fitted values are transformed into predicted probabilities.
- **absorb** (optional): allows for one dimension of fixed effects to be absorbed in the second-step estimation (in this case, the program relies on Stata's `areg`).
- **standardize** (optional): standardizes the coefficient on the generated regressor in the second-step estimation to allow for a direct comparison with the case in which the first-step dependent variable would be included directly as an independent variable in the second-step estimation. The standardization follows from equation (3) of [Hou and Loh (2016)](https://doi.org/10.1016/j.jfineco.2016.02.013).
- **fast** (optional): allows for much faster bootstrap execution using mata: the first-step regression estimates use the efficient computational algorithm of [MacKinnon (2021)](http://qed.econ.queensu.ca/pub/faculty/mackinnon/working-papers/qed_wp_1465.pdf), which speeds up the first-step estimations extremely; the second-step regression estimates cannot rely on this algorithm, but are also performed more efficiently given the use of mata functions.

---

Examples of output:

![image](https://user-images.githubusercontent.com/65561067/144615692-cc454959-ae79-4a40-9e5e-138c33a36eb0.png)

---

![image](https://user-images.githubusercontent.com/65561067/144615769-c49d1f8c-8fc7-44df-8ace-6f17c126e318.png)

---

![image](https://user-images.githubusercontent.com/65561067/144615851-fcaf8989-d714-4aef-b029-8d92c76ec601.png)

---

![image](https://user-images.githubusercontent.com/65561067/144615923-be200983-2d93-4ebe-a95f-fb962da48d7a.png)

---
