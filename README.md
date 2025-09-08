# SSMPG2025
Resources for [Software and Statistical Methods for Population Genetics 2025 (SSMPG 2025)](https://ssmpg2025.sciencesconf.org/) (Aussois, September 8-12 2025)


This public repository contains resources and instructions for the sixth edition of the Summer School on Software and Statistical Methods for Population Genetics (SSMPG 2025).
The SSMPG 2025 Summer School will provide a comprehensive overview of software and statistical methods for genotype-environment association analyses to detect genetic loci involved in local adaptation, as well as methods for predicting population maladaptation to future climate conditions.

##  Datasets

###  Snow slug
The dataset contains allele frequencies measured in 100 populations of the snow slug (*Limax nivalis*), an emblematic species in Olivier's computer that survived an abrupt climatic change. Sampled individuals were genotyped at L=2333 loci. Sixteen bioclimatic variables—comprising eight precipitation and eight temperature variables—were recorded both before and after the abrupt change in environmental conditions.

The objective of the practical session is to apply genotype-environment association methods to detect loci associated with fitness traits in the 100 snow slug populations. These loci will then be used to compute predictions of genetic offset for each population. The ultimate goal is to approximate as closely as possible the ground-truth measure of fitness loss observed in the populations.


###  Woolly marmot

The woolly marmot is an emblematic rodent species that lives in Thibaut's computer and recently disappeared from three directories due to poor coding decisions. Thibaut would like to repopulate these three sites but wants to ensure that the reintroduced individuals are optimally adapted to local environmental conditions.

The goal of the practical session is to select source populations that minimize maladaptation at the three reintroduction sites using genetic offset measures. The data consist of a genotype matrix for 610 individuals from 61 populations, genotyped at L=1000 loci (in 0/1/2 genotype format). Ten environmental variables have been measured for each source and reintroduction directory.

The associated datasets can be found in the Woolly marmot folder.



## Create teams, collaborate  and submit a paper

During the data analysis sessions, participants are encouraged to form teams. Ideally, each team should consist of 4–5 participants. Teams will collaboratively present a synthesis of their analysis in a public session on the final day and submit three files to the organizers for each data analysis challenge. The submitted files should include:

   * A list of candidate loci detected using their preferred GEA method, or a combination of methods, for both challenges

   * A list of genomic offset values obtained using their preferred offset method or a combination of methods for each snow slug population

   * A list of optimal source populations for the three Woolly marmot reintroduction sites

   * A short README file ("Materials and Methods") explaining all the methodological choices made during the analysis

Each team will be asked to present 2–3 slides for each dataset in the public session on the final day.

## 3. Evaluation

Thibaut and Olivier will reveal the truth about woolly marmots and snow slug, and the limits of inference from the genotypes and environmental data  Don't be worried. Everyone wins! 

##  4. Install software

### Install R and Rstudio
To participate in the practical sessions, bring your own laptop and install [R](https://cran.r-project.org/) and [RStudio](https://www.rstudio.com/), an integrated development environment (IDE) for R.

### Install R packages (LEA, gradientForest, vegan, qvalue)
To install R packages for the data analyses, copy and paste the following pieces of code in the R session

```r
#Install R packages for SSMPG 2025


#Package LEA (latest version) 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("LEA")

#Package gradientForest
install.packages("gradientForest", repos="http://R-Forge.R-project.org")

#Package vegan for RDA
install.packages("vegan")

#Package qvalue for controlling FDRs
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")

```
### Install BAYPASS

Download the archive for the latest stable version (3.1) from [here](https://forge.inrae.fr/mathieu.gautier/baypass_public) or directly via the following command run on a terminal:
```
wget https://forge.inrae.fr/mathieu.gautier/baypass_public/-/archive/v3.1/baypass_public-v3.1.tar.gz
```
Extract the archive, *e.g.*, from a terminal:
```
tar -zxvf baypass_public-v3.1.tar.gz
```
*BayPass* is coded in modern Fortran and can therefore be compiled for any system supporting a Fortran 90 (or higher) compiler. More specifically, the provided Makefile  is designed to work with:

- the free compiler *gfortran* (version >=10.0 recommended; if not already installed in your system, binaries may be available [here](https://gcc.gnu.org/wiki/GFortranBinaries) and are easy to install for most Windows, Mac and Linux OS versions) or;
- the commercial Intel® Fortran compilers *ifx* and *ifort* which are available at no cost for most platforms as
part of the [Intel® oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html).

Note also that *BayPass* relies on [*OpenMP*](https://www.openmp.org/wp/) to implement multi-threading (usually pre-installed on most Linux OS) on computer systems that have multiple CPUs or CPUs with multiple cores. 

With *gfortran*, *ifx* and *ifort*, compiling the source can simply be done with the following `make` instructions run within the *sources* sub-directory:

-  with *gfortran* (the command then produces an executable called `g_baypass`):
``` 
make clean all FC=gfortran
```

- with *ifx* (the command then produces an executable called `ifx_baypass`):
``` 
make clean all FC=ifx
```

- with *ifort* (the command then produces an executable called `i_baypass`):
``` 
make clean all FC=ifort
```
> Note: Under Linux (or MacOS), before the first use, make sure to give appropriate execution rights to the program. For instance you may run:
>```chmod +x baypass```

See the [manual](https://forge.inrae.fr/mathieu.gautier/baypass_public/-/blob/master/manual/BayPass_manual.pdf) for more details on how to compile and the program in general.



