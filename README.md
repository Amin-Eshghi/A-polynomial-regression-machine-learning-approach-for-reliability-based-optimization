# A-polynomial-regression-machine-learning-approach-for-reliability-based-optimization
This is a compact code for reliability analysis under uncertainty using a Polynomial Regression Machine Learning approach. The code implements a stochastic response surface method (SRSM) which quantifies the uncertainty in a performance function for the purpose of reliability-based design optimization (RBDO).
The method first uses the Polynomial Regression to approximate the limit state functions (LSF). Then it adds some weights and implements Moving Least Square method to locate the new design point closer to the estimated LSF. The two-step process results in a stochastic response surface, with which we can apply the Monte Carlo simulation (MCS) to obtain the full probabilistic characteristics (e.g., statistical moments, reliability, PDF and quantile) of the performance function. One mathematical example is solved for demonstration purpose.

FUNCTION DEFINITION

ML_RS(): main function
ML_sampling(): identify locations of polynomial regression samples
FindResponse(): define performance functions

VARIABLE DEFINITION

u: mean vector of random variables; s: standard deviation vector
cid: constraint number; rt: target reliability
nv: number of random variables; ns: number of MCS samples
uniComp: polynomial component function values
Response_RS: random response values generated using SRSM
GT_RS: quantile of response at target reliability (CDF) level (used as reliability constraint in RBDO)
Sen_GT_RS: sensivities of quantile w.r.t. means of design variables
To check out more details about this method and some numerical examples you can check out:

Eshghi, Amin Toghi, and Soobum Lee. "Adaptive improved response surface method for reliability-based design optimization." Engineering Optimization (2019): 1-19.
