# `gblc`: Boosted Lee-Carter

`gblc` provides a novel methodology for multi-population mortality forecasting proposed by Li Li, Han Li, and Anastasios Panagiotelis. They uses well-known stochastic mortality models as weak learners in gradient boosting rather than trees, and includes a penalty that shrinks the forecasts of mortality in adjacent age groups and nearby geographical regions closer together.


## Data
The proposed method demonstrates superior forecasting performance based on US male mortality data from 1969 to 2019.

## Usage
- For point forecasts, use `gb_forecast` and `gb_forecast_state` functions for national and state levels respectively.
- For prediction intervals, use `gb_forecast_pi` and `gb_forecast_pi_ss` functions for national and state levels respectively.
- The input mortality data should be organized as a T*n matrix. T denotes years and n denotes ages (for each state).
