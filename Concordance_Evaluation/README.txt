ConcordanceMetrics.m

Script: MATLAB R2024B release V 1.2
Author: L.I.G.J. july 2026.
Input= Change the values in the Contingency counts section


Reproduce Cohen's kappa and Gwet's AC1 for the provided 2x2 table
and verify with the Bootstrap method for small sample sizes.

Output= the results are displayed in the command Windows
	Observed agreement
	Cohen's kappa
	Approximate SE for kappa
	Kappa (z,p)
	Kappa CI(95%) approximation
	
	Gwet's AC1

	Kappa CI(95%) by Bootstrap
	AC1   CI(95%) by Bootstrap

References:
[1] Cohen, J. (1960). A coefficient of agreement for nominal scales. Educational and psychological measurement, 20(1), 37-46.
[2] McHugh, M. L. (2012). Interrater reliability: the kappa statistic. Biochemia medica, 22(3), 276-282.
[3] Gwet, K. L. (2008). Computing inter-rater reliability and its variance in the presence of high agreement. British Journal of Mathematical and Statistical Psychology, 61(1), 29-48.
[4] Sim, J., & Wright, C. C. (2005). The kappa statistic in reliability studies: use, interpretation, and sample size requirements. Physical therapy, 85(3), 257-268.
