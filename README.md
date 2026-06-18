# Atlantic Menhaden Mortality

Code and data supporting:

Ault, J.S. and Luo, J. (2026)

Investigation of Atlantic Menhaden Mortality Rates.  
The original code was from Liljestrand et al (2019)

## Repository Contents

Input data:
  Coston_ME85params.dat
  NMFS_ME70params.dat

ADMB model code:
  Coston_ME85params.tpl
  Coston_ME70params.tpl


## Requirements

ADMB 13.2+
R 4.5+

## Reproducing Results

1. Compile ADMB model

admb Coston_ME85params

2. Run model

   Coston_ME85params

3. Run MCMC

 Coston_ME85param -mcmc 4000000 -mcsave 1000
 Coston_ME85param -mceval

4. Run R scripts in scripts/R

## Contact

Jiangang Luo
University of Miami