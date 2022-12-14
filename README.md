# DNS-merger-rate-2022
Modifications to PsrPopPy2 used for the MMGPS-L PSR J1208-5936 Bernadich et al. 2023 paper.

## Description
This is an expansion to the code used in K. Gruntal et al. 2021, where a new pulse model for J1906+0746 is implemented. In thurn, that code is an expansion from https://github.com/NihanPol/2018-DNS-merger-rate

In this implementation, all pulsars are callable from the single ```find_alpha_general.py``` scripts, and for the MMGPS-L survey a new degradation factor has been implemented to ```dosurvey.py``` to account for the interferometric nature of MeerKAT. This is explained in section 6 of Bernadich et al. 2023

## Usage
Just substitute the files from your PsrPoPy2 installation (https://github.com/devanshkv/PsrPopPy2) with the ones provided here, and add MMGPS-L to the surveys.
