# Simulation of Markov Jump Processes

Lab work for the course *Stochastic Models of Dynamical Systems*, Applied Mathematics (4th year), INSA Rennes, 2026.
Team project with Kelian Esnoult.

## Goal

Markov jump processes model systems that move between states at random times in continuous time: customers arriving
in a queue, machines breaking down, a borrower's credit rating changing. This project builds simulation and computation
tools for these processes and checks them against theory.

## What we did

1. **Poisson process**: simulation of event times over an interval, and check that the waiting times between events
   follow the expected exponential distribution.
2. **M/M/K queues** (customers served by K servers): simulation of the number of customers over time in the stable,
   critical and explosive regimes. In the stable regime, the long-run distribution of the number of customers matches
   the theoretical one, confirmed by a chi-squared goodness-of-fit test.
3. **Uniformization**: computation of the probability of being in each state at any time *t* (the matrix exponential
   of the generator), with a bound on the numerical error. We use it to detect the time from which the process has
   reached its long-run regime.
4. **Conditioned paths (Markov bridges)**: simulation of trajectories forced to start and end in given states, or to pass
   through observed states at given times, as needed when a process is only observed at discrete dates.

## Link with finance

The same tools appear in quantitative finance: Poisson processes drive jump models of asset prices and default-intensity
models in credit risk, and the matrix exponential of a generator is used to derive credit rating migration probabilities
over any horizon.

## Files

- `Rapport_TP1_JMP_Yann_Robcis_Kelian_Esnoult.pdf`: full report with proofs and results (in French)
- `Code_TP1_JMP_Yann_Robcis_Kelian_Esnoult.r`: R code
