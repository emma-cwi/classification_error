######### BASIC FUNCTIONS
code_biasCorrection.R
Bias correction methods: misclassificatino method, reclassification method, Ratio-to-TP method

code_fieller.R
Basic function to calculate confidence intervals from Fieller's theorem, i.e., CI for ratio A/B of correlated variables A and B.
Functions to apply Fieller's interval to estimates of 1) class size nx., 2) detailed error composition nxy, and 3) reproduction of prior work in Shieh 2009 and Buonaccorsi 2010. Applications of Fieller's theorem for estimates 1) and 2) use Sample-to-Sample method to estimate the variance of variables A and B.


######### TEST FUNCTIONS
test_allMethods.R
For set of tests as in Beauxis-Aussalet 2017

code_testBiasCorrection.R
Tests for bias correction methods.

code_testFieller.R
Test for Sample-to-Sample methods, applied to derive confidence intervals using Fieller's theorem.

code_testMaxDeterminant.R
Tests for Maximum Determinant method.


######### DEFINITION OF COMMON VARIABLES
test_n0., test_n1.
True class size in test set (for class 0 and 1, respectively).

target_n.0, target_n.1
"Imperfect" class size in target set, as outputed by classifier.

misclass_n0., misclass_n1.
"Corrected" class size in target set, resulting from misclassification method.

t01, t10
Error rates as measured in test set, respectively, t01 = 1 - TP Rate and t10 = FP Rate.

nTest, nTarget
Number of test and target sets (respectively) sampled for testing the methods.

z
Intended coverage of confidence intervals, i.e., z=1 for 68% CI for , z=1.96 for 95% CI.
