# STI Screening among Men who have Sex with Men (MSM) in the United States

This repository holds the source to code to reproduce the analysis featured in our STI and HIV transmission model among men who have sex with men in the United States. This study investigated the possible scale-up of the CDC's screening recommendations for MSM in the United States. Different scale-up strategies included screening coverage, screening frequency, and definition of higher-risk sexual behavior, including the impact of these strategies on STI incidence among MSM.

## Citation

> Weiss KM, Jones JS, Anderson EJ, Gift TL, Chesson H, Bernstein K, Workowski K, Tuite A, Rosenberg ES, Sullivan PS, Jenness SM. Optimizing Coverage vs Frequency for Sexually Transmitted Infection Screening of Men Who Have Sex With Men. _Open Forum Infect Dis._ 2019; 6(10): ofz405.

<img src="https://github.com/EpiModel/stitestguidelines/raw/master/analysis/Fig1.png">

## Abstract
Incidence of bacterial sexually transmitted infections (STIs) in men who have sex with men (MSM) has increased substantially despite availability of effective antibiotics. The US CDC recommends annual screening for all sexually active (SA) MSM and more frequent screening for high-risk (HR) MSM. The population-level benefits of improved coverage versus increased frequency of STI screening among SA versus HR MSM are unknown. We used a network transmission model of gonorrhea (NG) and chlamydia (CT) among MSM to simulate the implementation of STI screening across different scenarios, starting with the CDC guidelines at current coverage levels. Counterfactual model scenarios varied screening coverage and frequency for SA MSM compared to HR MSM, who had multiple recent partners. We estimated infections averted and the number needed to screen to prevent one new infection. Compared with current recommendations, increasing the frequency of screening to biannual for all SA MSM and adding 5%some HR screening coverage could avert 72% of NG and 78% of CT infections over 10 years. Biannual screening of 30% of HR MSM, at empirical coverage levels for annual SA screening, could avert 76% of NG and 84% of CT infections. Other scenarios, including higher coverage among SA MSM and increasing frequency for HR MSM, averted fewer infections but did so at a lower number needed to screen. The optimal screening scenarios in this model to reduce STI incidence among MSM included more frequent screening for all sexually active MSM and higher coverage of screening for HR men with multiple partners.

<img src="https://github.com/EpiModel/stitestguidelines/raw/master/analysis/Fig2.png">

## Model Code

These models are written and executed in the R statistical software language. To run these files, it is necessary to first install our epidemic modeling software, [EpiModel](http://epimodel.org/), and our extension package specifically for modeling HIV and STI transmission dynamics among MSM, [EpiModelHIV](http://github.com/statnet/EpiModelHIV).

In R:
```
install.packages("EpiModel", dep = TRUE)

# install remotes if necessary, install.packages("remotes")
remotes::install_github("statnet/tergmLite")
remotes::install_github("statnet/EpiModelHPC")
remotes::install_github("statnet/EpiModelHIV", ref = "ept")
