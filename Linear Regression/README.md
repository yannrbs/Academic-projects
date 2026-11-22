# Forecasting US Air Traffic from Macroeconomic Data

Linear regression project, Applied Mathematics (3rd year), INSA Rennes, 2026.
Team project with Kelian Esnoult.

## Goal

Explain the annual number of passengers carried by US airlines with US macroeconomic indicators,
then forecast traffic over 2026–2034.

**Data:** 29 annual observations (2006–2034). The traffic is known up to 2025; for 2026–2034 only the
economic variables are available.

| Variable | Description |
|----------|-------------|
| TotalPassengers_US | Passengers carried per year by US airlines (variable to explain) |
| CPI | Consumer Price Index |
| Conso_gouv | Government consumption |
| Conso_privee | Private consumption |
| Exports / Imports | US exports and imports |
| GDP | Gross Domestic Product |
| Investissements | Fixed investment |

## Method

1. **Log transformation** of all variables, to stabilise the variance and make the relationships more linear.
2. **Backward stepwise selection**: non-significant variables are removed one at a time, and each reduction is
   checked with a Fisher test against the full model.
3. **Residual diagnostics**: residual plots, normality (Q-Q plot), Durbin-Watson test for autocorrelation,
   leverage and Cook's distance to detect influential points.
4. **Three iterations**: 2009 (financial crisis), 2020 and 2021 (Covid and recovery) are removed one after the other,
   as exceptional shocks that do not reflect the usual behaviour of air traffic.
5. **Forecasting**: predictions for 2026–2034 with 95% prediction intervals, converted back from logs.

## Results

The final model keeps four variables, all significant, with an adjusted R-squared of 0.975:

```
log(Passengers) = -6.403 - 1.547 log(CPI) + 2.002 log(Conso_gouv) + 0.695 log(Exports) + 0.868 log(GDP)
```

- Residuals are centred, close to normal, with stable variance and no autocorrelation (Durbin-Watson p-value 0.587).
- The signs make economic sense: higher prices reduce traffic, while growth in GDP, exports and public spending increases it.
- The model forecasts a steady rise in traffic over 2026–2034.

**Limitations:** the model is fitted on only 17 observations after removing the crisis years, the explanatory variables
for 2026–2034 are themselves forecasts, and the model cannot anticipate future shocks. The prediction interval widens
with the horizon, which reflects this uncertainty.

## Files

- `Robcis_Esnoult_Rapport_Projet_Regression_Lineaire.pdf`: full report (in French)
- `projet.r`: R code
