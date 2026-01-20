([Français](#beast_pype-))

# **BEAST_pype**: An automated pipeline for high throughput phylodynamic analyses using BEAST 2.

BEAST_pype is a pipeline with the aim of automating and parallelizing many of the steps involved in phylodynamics using 
[BEAST 2](https://www.beast2.org/). The [workflows](https://github.com/m-d-grunnill/BEAST_pype/wiki) 
expedite phylodynamic analyses using [BEAST 2](https://www.beast2.org/) via a series of python-based Jupyter Notebooks. 

BEAST_pype was concieved with the aim to speed up and automate the use of BEAST-2 analyses for routine public health use at the Public Health Agency of Canada.  
Activities include running routine analyses on continuously circulating viruses, such as SARS-CoV-2 or Influenza, to extract epidemiological parameters of interest for surviellance. Also, for accelerating research or outbreak investigations that require several experimental runs for optimization and more rapid results generation.  

Features of BEAST_pype include:  
* Command-line launchable, and Yaml file controlled pipeline.
* Improved XML generation and re-use, e.g. using user-provided BEAST-2 xmls as a template for generating a new BEAST-2 xml
   from new sequences with the associated metadata.
* Launching several parallelized runs and a GUI for selecting converged MCMC chains.
* Automated generation of reports analysing BEAST-2 runs (diagnostic and result plots and statistics).
* Ability to create and include initial trees for faster runs.
* Incorporating downsampling techniques.
  
### For instructions, getting started, and more details, please see [BEAST_pype's wiki](https://github.com/m-d-grunnill/BEAST_pype/wiki).

## IMPORTANT NOTES
* For ease of distribution reasons beast_pype uses the version of BEAST-2 that is currently available via conda, specifically [bioconda](https://anaconda.org/bioconda/beast2), 2.6 as 2025-Jun-25. Template BEAST-2 XMLs from other versions of may not work. BEAST 2.7.7 is available on the conda channel [millerjeremya](https://anaconda.org/millerjeremya/beast2). However, this has been tested on a Linux OS (2025-06-04) and the command line arguments did not work.

## Installation instructions:

See the [wiki's installation instructions section](https://github.com/m-d-grunnill/BEAST_pype/wiki/Installation-Instructions). 
   
## Guides on Running BEAST_pype Workflows

The documentation for running BEAST_pype's workflows is in the repo's [wiki](https://github.com/m-d-grunnill/BEAST_pype/wiki).

## How to Contribute

See [CONTRIBUTING.md](CONTRIBUTING.md)

## Legal
           

> Copyright (c) His Majesty the King in Right of Canada, as represented by the 
> Minister of Health, 2025.

Unless otherwise noted, the source code of this project is covered under Crown Copyright, Government of Canada, and is distributed under the [GPL-2 license](LICENSE).

The Canada wordmark and related graphics associated with this distribution are protected under trademark law and copyright law. No permission is granted to use them outside the parameters of the Government of Canada's corporate identity program. For more information, see [Federal identity requirements](https://www.canada.ca/en/treasury-board-secretariat/topics/government-communications/federal-identity-requirements.html).

## Contributors

Martin Grunnill [martin.grunnill@phac-aspc.gc.ca](mailto:martin.grunnill@phac-aspc.gc.ca)  
Carmen Lia Murall [carmen.lia.murall@phac-aspc.gc.ca](mailto:carmen.lia.murall@phac-aspc.gc.ca)   
Rachelle Di Tullio   
Kodjovi Mlaga  

## Acknowledgments

### Useful Software

A forerunner and inspiration of this work was [CoV flow](https://gitlab.in2p3.fr/ete/CoV-flow).

BEAST_pype would not be possible without:
* [BEAST 2](https://www.beast2.org/) also [Bouckaert R. *et al.* (2019)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006650)
* The initial work on [BEAST2xml](https://github.com/acorg/beast2-xml/tree/master) by Terry Jones.

### Feedback During Development

George Long  
Venketa Duvvuri  
Louis Du Plessis  
Samuel Alizon  
Remco Bouckaert
Gonche Danesh

______________________

# BEAST_pype: ...

- Quel est ce projet?
- Comment ça marche?
- Qui utilisera ce projet?
- Quel est le but de ce projet?

## Comment contribuer

Voir [CONTRIBUTING.md](CONTRIBUTING.md)

## Licence

Sauf indication contraire, le code source de ce projet est protégé par le droit d'auteur de la Couronne du gouvernement du Canada et distribué sous la [licence GPL-2](LICENSE).

Le mot-symbole « Canada » et les éléments graphiques connexes liés à cette distribution sont protégés en vertu des lois portant sur les marques de commerce et le droit d'auteur. Aucune autorisation n'est accordée pour leur utilisation à l'extérieur des paramètres du programme de coordination de l'image de marque du gouvernement du Canada. Pour obtenir davantage de renseignements à ce sujet, veuillez consulter les [Exigences pour l'image de marque](https://www.canada.ca/fr/secretariat-conseil-tresor/sujets/communications-gouvernementales/exigences-image-marque.html).




