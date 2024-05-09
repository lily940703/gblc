## For reproducibility purpose

#### Environment requirements

R >= 4.2.2 with the following required packages:

- `demography` (1.22)
- `lightgbm` (3.3.5)
- `MCS` (0.1.3)
- `ggplot2` (3.4.2)
- `reshape2` (1.4.4)
- `splm` (1.6-2)
- `parallel` (4.2.2)
- `doParallel` (1.0.17)
- `foreach` (1.5.2)

### Usage

We organize the code that replicates all the results in the paper：

- `forecast_national.R` replicates Tables 4.1, 4.2, B.1, C.1 and Figures 4.1, 4.2, 4.3.
- `forecast_state.R` replicates Tables 4.3, 4.4, 4.5, A.2, B.2, C.2 and Figures 4.4, 4.5, 4.6, 4.7.

Because processing state-level data is time-consuming, we provide intermediate results (/data) that can be loaded if needed.

When running the code, please change the path "D:/gblc-main/" to the real path on your computer.

