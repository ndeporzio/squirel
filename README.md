# SquiRel - Squishable Distribution Relics
**Authors:** Nick Deporzio, David Imig, Jessie Shelton

We provide our modifications to CLASS and MontePython to reproduce our results and plots of the SquiRel paper, and detail our custom power-law sampling algorithm for handling logarithmic distributions. 

## Instructions
- Set up SquiRelCLASS following the standard [CLASS](https://github.com/lesgourg/class_public) installation
- Run the `main.py` plotting script and specify which plots are desired when prompted
- Depending on the desired plot, the code may run a number of CLASS instances
- To recreate our contour plots, a working [GetDist](https://getdist.readthedocs.io/en/latest/) is necessary, and one must download the MCMC chains from [our DropBox](https://www.dropbox.com/scl/fo/pljewqsxyg1hzwqugb7he/AGTLiscyh-99CvWGXB5EXw4?rlkey=arfsaa2620y4mg8gehgniz26h&st=hf8iie9h&dl=0)

## Power Law Sampling Scheme
Our custom sampling scheme is implemented in `quadrature.c` of SquiRelCLASS. It captures specific details about the shape of an unwieldy distribution near its maximum while also adequately sampling its potentially significant momentum tails, all while using many fewer momentum samples than CLASS' default sampling strategies. 

In detail, given a distribution of the form
```math
f_\mathrm{LN}(\xi\mid\sigma_\mathrm{LN})=\frac{ 1 }{ \xi \sigma_\mathrm{LN} \sqrt{2\pi } } \exp\left( - \frac{ \ln \left( \xi \right)^2}{ 2 \sigma_\mathrm{LN}^2 } \right),
```
that we wish to integrate over $\xi$, linearly partitioning the range between $\xi_{\rm min}$ and $\xi_{\rm max}$ is an inefficient strategy. We instead re-scale the integration variable according to a power-law of the form 
```math
\xi \rightarrow \xi' \equiv \xi^\alpha,
```
which can accurately sample the distribution using fewer sampling points. For each choice of distribution parameter $\sigma$, optimal values for $\alpha$, $\xi_{\rm min}$, $\xi_{\rm max}$, and the number of sampling points $N$ must be identified. The procedure we follow to numerically identify these optimal values is:
1. Establish an accuracy threshold $\delta_{\rm acc} = 10^{-6}$,
2. At each $\sigma_{\rm LN}$, find $\xi_{\rm min}$ and $\xi_{\rm max}$ such that $f_{\rm LN}(\xi_{\rm min}\mid\sigma_\mathrm{LN}) = f_{\rm LN}(\xi_{\rm max}\mid\sigma_\mathrm{LN})$ and $\delta_{\rm analytic}<\delta_{\rm acc}$ where $\delta_{\rm analytic} \equiv |(\int_{\xi_{\rm min}}^{\xi_{\rm max}} f_{\rm LN}(\xi\mid\sigma_\mathrm{LN})d\xi - \int_0^\infty f_{\rm LN}(\xi\mid\sigma_\mathrm{LN})d\xi) / \int_0^\infty f_{\rm LN}(\xi\mid\sigma_\mathrm{LN})d\xi |$.
3.  For each $\alpha \in [0.1, 0.2, 0.3, ..., 100]$ find $N_\alpha$,  the smallest number of $\xi_i$ sampled from $f_\mathrm{LN}(\xi\mid\sigma_\mathrm{LN})$ such that $\delta_{\rm numeric}<\delta_{\rm acc}$ where $\delta_{\rm numeric}\equiv|\left(\sum_{\xi_{\rm min}}^{\xi_{\rm max}} f_{\rm LN}(\xi\mid\sigma_\mathrm{LN})\Delta\xi - \int_0^\infty f_{\rm LN}(\xi\mid\sigma_\mathrm{LN})d\xi\right) / \int_0^\infty f_{\rm LN}(\xi\mid\sigma_\mathrm{LN})d\xi|$
5.  identify that $\alpha$ which satisfies $\delta_{\rm numeric}<\delta_{\rm acc}$ with the smallest $N_\alpha$ as the optimal choice of $\alpha$ for the given $\sigma_{\rm LN}$,
6.  fit the $\alpha(\sigma_{\rm LN})$ points to a third degree polynomial yielding $\alpha_{\rm fit}(\sigma_{\rm LN})$,
7.  For each $\alpha_{\rm fit}(\sigma_{\rm LN})$ we re-compute the necessary number of $\xi'$ points $N_{\alpha, \textrm{fit}}$ needed to achieve $\delta_{\rm numeric}<\delta_{\rm acc}$,
8. we identify as $N$ the largest $N_{\alpha, \textrm{fit}}$ across all values of $\alpha_{\rm fit}(\sigma_{\rm LN})$,
9. finally, for each choice of $\sigma_{\rm LN}$, we numerically sample $f_\mathrm{LN}(\xi\mid\sigma_\mathrm{LN})$ according to the $\xi \rightarrow \xi' \equiv \xi^\alpha$ rescaling using $\xi_{\rm min}$, $\xi_{\rm max}$, $\alpha_{\rm fit}(\sigma_\mathrm{LN})$ and $N$.  

To validate our custom sampling scheme, we confirm that maximum relative differences between the CMB spectra obtained by sampling LN distributed LiMRs using our scheme and the CLASS default scheme remain less than 0.1 percent for $\sigma\in\left[0.04,1.5\right]$. 

## Contact
Please email deporzio@bu.edu, dimig2@illinois.edu, or jshelton137@gmail.com for any questions.

## Publications
If you make use of this code in your publication, please cite our paper along with CLASS and (possibly) MontePython and GetDist. 

```
% CLASS:
@article{Blas:2011rf,
    author = "Blas, Diego and Lesgourgues, Julien and Tram, Thomas",
    title = "{The Cosmic Linear Anisotropy Solving System (CLASS) II: Approximation schemes}",
    eprint = "1104.2933",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    reportNumber = "CERN-PH-TH-2011-082, LAPTH-010-11",
    doi = "10.1088/1475-7516/2011/07/034",
    journal = "JCAP",
    volume = "07",
    pages = "034",
    year = "2011"
}
% MontePython:
@article{Brinckmann:2018cvx,
      author         = "Brinckmann, Thejs and Lesgourgues, Julien",
      title          = "{MontePython 3: boosted MCMC sampler and other features}",
      year           = "2018",
      eprint         = "1804.07261",
      archivePrefix  = "arXiv",
      primaryClass   = "astro-ph.CO",
      SLACcitation   = "%%CITATION = ARXIV:1804.07261;%%"
}
@article{Audren:2012wb,
      author         = "Audren, Benjamin and Lesgourgues, Julien and Benabed,
                        Karim and Prunet, Simon",
      title          = "{Conservative Constraints on Early Cosmology: an
                        illustration of the Monte Python cosmological parameter
                        inference code}",
      journal        = "JCAP",
      volume         = "1302",
      pages          = "001",
      doi            = "10.1088/1475-7516/2013/02/001",
      year           = "2013",
      eprint         = "1210.7183",
      archivePrefix  = "arXiv",
      primaryClass   = "astro-ph.CO",
      reportNumber   = "CERN-PH-TH-2012-290, LAPTH-048-12",
      SLACcitation   = "%%CITATION = ARXIV:1210.7183;%%",
}
% GetDist
@article{Lewis:2019xzd,
 author         = "Lewis, Antony",
 title          = "{GetDist: a Python package for analysing Monte Carlo
                   samples}",
 year           = "2019",
 eprint         = "1910.13970",
 archivePrefix  = "arXiv",
 primaryClass   = "astro-ph.IM",
 SLACcitation   = "%%CITATION = ARXIV:1910.13970;%%",
 url            = "https://getdist.readthedocs.io"
}
```
