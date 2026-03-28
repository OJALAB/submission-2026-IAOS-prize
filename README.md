# Simulation Code: Submission for the 2026 IAOS Prize for Young Statisticians (Adam Struzik)

This repository contains the simulation code for the paper *Correcting for Misclassification and Missing Data in Job Vacancy Estimation: A Capture-Recapture Approach*
submitted for the 2026 IAOS Prize for Young Statisticians by Adam Struzik.

---

## Project Structure

```
submission-2026-IAOS-prize/
├── code/                       
    ├── functions.R             Simulation core functions
    └── simulation.R            Main simulation script
├── data-raw/                   
    └── data.RData              Data from ePraca
├── results/
    ├── results_total.RData     Estimates of the total number of vacancies
    └── results_by_code.RData   Estimates of the number of vacancies for different occupations
```

## Usage

To replicate the results, run `code/simulation.R`.