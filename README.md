# SquiRel - Squishable Distribution Relics
**Authors:** Nick Deporzio, David Imig, Jessie Shelton

We provide our modifications to CLASS and MontePython to reproduce our results and plots of the SquiRel paper. 

## Instructions
- Set up SquiRelCLASS following the standard [CLASS](https://github.com/lesgourg/class_public) installation
- Run the `main.py` script and specify which plots are desired when prompted
- Depending on the desired plot, the code may run a number of CLASS instances
- If contour plots are desired, a working [GetDist](https://getdist.readthedocs.io/en/latest/) is necessary, and one must download our MCMC chains from [our DropBox] (https://www.dropbox.com/scl/fo/pljewqsxyg1hzwqugb7he/AGTLiscyh-99CvWGXB5EXw4?rlkey=arfsaa2620y4mg8gehgniz26h&st=hf8iie9h&dl=0)

## Contact
Please email deporzio@bu.edu, dimig2@illinois.edu, or jshelton137@gmail.com for any questions.

## Publications
If you make use of this code in your publication, please cite our paper along with MontePython, GetDist and CLASS. 

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
