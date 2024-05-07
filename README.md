# `gblc`: Boosted Lee-Carter

`gblc` provides a novel methodology for multi-population mortality forecasting proposed by Li Li, Han Li, and Anastasios Panagiotelis. They uses well-known stochastic mortality models as weak learners in gradient boosting rather than trees, and includes a penalty that shrinks the forecasts of mortality in adjacent age groups and nearby geographical regions closer together.


## Data
The proposed method demonstrates superior forecasting performance based on US male mortality data from 1969 to 2019.

## Usage
- For point forecasts, use `gb_forecast` and `gb_forecast_state` functions for national and state levels respectively.
- For prediction intervals, use `gb_forecast_pi` and `gb_forecast_pi_ss` functions for national and state levels respectively.
- The input mortality data should be organized as a T*n matrix. T denotes years and n denotes ages (for each state).

## Quick Demo for national point forecasts

### Data preparation

```{r}
data = read.table(".../Mortality_national_male.txt",head=T)
data_y = data[data[,"Year"]%in%c(1969:2019),]
data_pop = read.table(".../Population_national.txt",head=T)
data_pop_y = data_pop[data_pop[,"Year"]%in%c(1969:2019),]

# T*n matrix
matrix_motality_male = c()
for (i in 1969:2019) {
  row =  data_y[(data_y[,"Year"]==i & data_y[,"Age"]%in%c(0:84)),"Male"]
  row_pop = data_pop_y[(data_pop_y[,"Year"]==i & data_y[,"Age"]%in%c(0:84)),"Male"]
  motality_comb = data_y[(data_y[,"Year"]==i & !data_y[,"Age"]%in%c(0:84)),"Male"]
  pop_comb = data_pop_y[(data_pop_y[,"Year"]==i & !data_pop_y[,"Age"]%in%c(0:84)),"Male"]
  age_comb = sum(motality_comb*pop_comb)/sum(pop_comb)
  row_comb = c(row, age_comb)
  matrix_motality_male = rbind(matrix_motality_male, row_comb)
}
colnames(matrix_motality_male) = c(0:84,"85+")
rownames(matrix_motality_male) = c(1969:2019)

```

### W matrix that captures the age structure

```{r}
matrix_W=diag(c(0,1,rep(2,82),1,0))
for (i in 2:85) {
  for (j in 2:85) {
    if(abs(i-j)==1){
      matrix_W[i,j]=-1
    }
  }
}
```

### Forecast 2010 to 2019
```{r}
  matrix_motality_train = matrix_motality_male[c(1:41),]
  matrix_motality_test = matrix_motality_male[c(42:51),]
  # GBLC: set lambda=0;
  # GBLC-age: set lambda>0;
  fore_GB_lca = gb_forecast(matrix_input=matrix_motality_train, 
                            matrix_W = matrix_W, 
                            lambda=0, h=10,
                            max.interation=100,scale=F,threshold=1e-06)
  mase_gblc = c()
  for (h in 1:10) {
    mase = get_mase(fore_GB_lca$final_forecast[c(1:h),], 
                    matrix(matrix_motality_test[c(1:h),],nrow=h), 
                    matrix_motality_train) 
    mase_gblc = c(mase_gblc,mase)
  }
  print("Out-of-sample mortality forecasts for years 2010-2019")
  fore_GB_lca$final_forecast
  print("MASE for horizons 1-10")
  mase_gblc

```

## For reproducibility purpose

#### Environment requirements

R >= 4.2.2 with the following required packages:

- `demography` (1.22)
- `lightgbm` (3.3.5)
- `MCS` (0.1.3)
- `ggplot2` (3.4.2)
- `reshape2` (1.4.4)


