To run this example, copy its content to a directory from which you can run seqEM.
>sequencer one.ctrl
It should provide a file identical to sample_result.txt

r = referent
v = variant

The individuals are ordered the same way as the input files in the control file. So:

[0] = one_0.out
[1] = one_1.out
[2] = one_2.out
[3] = one_3.out
[4] = one_4.out
[5] = one_5.out
[6] = one_6.out
[7] = one_7.out
[8] = one_8.out
[9] = one_9.out

The summary of each position/SNP is listed at the beginning:

0. 385856  iter= 24  error= 0.05329  p= 0.47001 hwe-chisq= -1.000000
P(rr)= 0.21950 p(rv)= 0.50101 p(vv)= 0.27949

0. 				the 0th position
385856			SNP id number
iter= 24			number of iterations for the EM-alogrithm to converge
error= 0.05329		estimated genotyping error. A high error usually indicates a problem.
p= 0.47001			referent allele frequency
hwe-chisq= -1.000000	Gives estimate of HWE equilibrium if at minimum 10 individuals are typed
P(rr)= 0.21950 p(rv)= 0.50101 p(vv)= 0.27949	estimated likelihood of genotypes

Indvidual results for [0]:
genotype[0]= rv props(rr,rv,vv)=  0.000  0.944  0.056    reads=     10 , variants=      7

Most likely genotype "rv".
Probability of genotypes rr, rv and vv are 0.000  0.944  0.056 respectively.
The read depth (= 100 and the number of variant reads are listed last.


