# Mathematical Modeling and Health-Economic Evaluation of Hepatitis B Vaccination in China

This repository contains the MATLAB source code and scripts for the dynamic compartmental modeling, MCMC parameter estimation, and health-economic evaluation (Cost-Benefit and DALYs) of the Hepatitis B virus (HBV) vaccination program in China.

## Repository Structure

- **`mcmc_function.m`, `mcmc_model.m`, `mcmc_ss.m`**: Scripts and functions for MCMC parameter estimation (DRAM algorithm) using epidemiological data.
- **`C1.m`, `C2.m`, `C3.m`**: Long-term population and infection projection scripts under different intervention scenarios (Monte Carlo simulations).
- **`HBV1.m`, `HBV2.m`, `HBV3.m`**: Fine-stage disease progression and mortality projection scripts.
- **`process_data.m` / `process_data_2.m`**: Data preprocessing and structural formatting scripts.
- **`d.m` / `d_2.m`**: Main analysis scripts for generating health-economic evaluation outputs (ICER, BCR, DALYs) and comparison tables.
- **`Sensitivity_analysis.m`**: Master script for automated sensitivity analyses across varying parameters.

## Requirements

- MATLAB (Tested on R2021a or later)
- MCMCstat toolbox (required for Bayesian inference scripts)
- Standard MATLAB Toolboxes (Statistics and Machine Learning Toolbox)

## Usage Instructions

1. Ensure all input datasets are placed in the `datasets/` directory.
2. Run the MCMC estimation and baseline scripts in `results/C1/`.
3. Run the projection models (`HBV1.m`, etc.) to generate fine-stage outcomes.
4. Execute `d.m` to reproduce the final tables and economic evaluation metrics.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
