# Simulation Code: Submission for the 2026 IAOS Prize for Young Statisticians (Adam Struzik)

This repository contains the simulation code for the paper *Estimation of Job Vacancy Statistics Under Misclassification and Missing Data Using Capture-Recapture*
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
The simulation is currently configured for parallel execution on 10 workers, so at least 10 CPU cores are recommended.