# Template BEAST 2 XMLs

This directory houses some [template-xml](Glossary-of-Terms.md#template-xml)s for use with BEAST_pype workflows.
Below are some notes on the specific templates.

## [BDSKY_serial_COVID-19_template.xml](BDSKY_serial_COVID-19_template.xml)

Birth Death Skyline Serial (BDSky serial) template for COVID-19.

Clock Model:
* Stict Clock gamma distribution Mean=4e-4 SD=3e-4 {'a': 1.7777777777777781, 'scale': 0.00022499999999999994} see [COVID_epi_params.ipynb](COVID_epi_params.ipynb).

Site Model:
* Gamma Catagory 4
* Shape 0.5 default prior
* Subs Model HKY 
* Kappa 2 default prior
* Frequencies default prior


BDSky serial Model priors:
* Rate becoming uninfectious Gamma {'a': 5.921111111111111, 'scale': 12.32876712328767, 'loc': 0} (see 5 days section of [New_COVID_inf_period.ipynb](New_COVID_inf_period.ipynb)), initial 36.25, boundaries [0, 365].
* Origin Uniform [0, 6] initial 6. This initial value is only so high as BEAST can generate an initial tree of height 2-5 and the origin initial value needs to be greater than that initial tree height.   
* $R_e$ Log Normal mean (not real) 1, s 0.5, initial 2.0, boundaries [0, 15]
* Sampling proportion Beta [1, 999], initial 0.01, boundaries [0, 1]

## [covid-19_coalescent_exponential.xml](covid-19_coalescent_exponential.xml)

Clock Model:
* Stict Clock gamma distribution Mean=4e-4 SD=3e-4 {'a': 1.7777777777777781, 'scale': 0.00022499999999999994} see [COVID_epi_params.ipynb](COVID_epi_params.ipynb).

Site Model:
* Gamma Catagory 4
* Shape 0.5 default prior [0.1 lower bound]
* Subs Model HKY 
* Kappa 2 default prior
* Frequencies default prior

Coalescent Exponential
ePopSize lognormal {'m': 1, 's' 5.0} initial = 3.
growth rate laplace distribution {mu: 0, scale: 100} initial = 3e-4.
 
