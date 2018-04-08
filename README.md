EE 511 Simulation Methods for Stochastic Systems Project #4: Investigations on Monte Carlo Methods
==============================================
[Pi-Estimation]
----------------
Generate n=100 samples of i.i.d 2-dimensional uniform random variables in the unit-square. Count how
many of these samples fall within the quarter unit-circle centered at the origin. This quarter circle
inscribes the unit square as shown below:

i.) Use these random samples to estimate the area of the inscribed quarter circle. Use this area estimate
to estimate the value of pi. Do k=50 runs of these pi-estimations. Plot the histogram of the k=50 piestimates
(each estimate based on n=100 2-D uniform samples).

ii.) Repeat the experiment with different numbers of uniform samples, n. Plot the sample variance of the
pi-estimates for these different values of n. Keep k=50 for all these runs. Comment on the sample
variance of your estimates.

[Monte Carlo Integration and Variance Reduction Strategies]
----------------------------------------------------------------
Use n=1000 random samples to obtain Monte Carlo estimates for the definite integrals:

(a) [1+sinh(2x)ln(x)]-1, x in [0.8,3]

(b) Exp[-x^4 – y^4], (x, y) in [-pi, pi]

Calculate the sample variance of the Monte Carlo estimates using a similar method as in problem 1.
Use the same number of random samples, n=1000, to obtain those Monte Carlo estimates. But this time
incorporate stratification and importance sampling in the Monte Carlo estimation procedures. Compare
the Monte Carlo estimates and their sample variances.

Discuss the quality of the Monte Carlo estimates from each method. Also discuss the strengths and
weaknesses of stratification and importance sampling in Monte Carlo estimation.

[Test your integral estimator on the following function with your own choice of n samples]
--------------------------------------------------------------------------------------------
𝑓(𝑥, 𝑦) = 20 + 𝑥^2 + 𝑦^2 − 10(𝑐𝑜𝑠[2𝜋 × 𝑥] + 𝑐𝑜𝑠[2𝜋 × 𝑦])

(x, y) in [-5, 5] for f(x,y)

[Turn in]:
-----------
- A summary of your experiments including histograms of your data, any relevant plots, and
relevant numerical results
- brief discussions of the results
- a print out of your code.
