Written by Chun-Lin Liu
E-mail: cl.liu@caltech.edu
Project website:
http://systems.caltech.edu/dsp/students/clliu/coarray.html
Last revised on July 13, 2017


1) Please execute "main_coarray_MUSIC_interpolation.m" in MATLAB. 
	cvx package (http://cvxr.com/cvx/) is required in Line 71.
	It is expected to obtain two figures similar to Figure 2 of [1].

2) weight_function.m returns the weight function of an array.

3) sample_covariance_to_difference_coarray.m [1,2] takes the sample covariance matrix and computes sample autocorrelation vectors or augmented covariance matrices on the difference coarray.

4) coarray_MUSIC.m [1] computes the coarray MUSIC spectrum or estimates the normalized DOAs.


References
[1] C.-L. Liu and P. P. Vaidyanathan, “Remarks on the Spatial Smoothing Step in Coarray MUSIC,” IEEE Signal Processing Letters, vol. 22, no. 9, pp. 1438-1442, Sep. 2015. 
[2] C.-L. Liu, P. P. Vaidyanathan and P. Pal, “Coprime Coarray Interpolation for DOA Estimation via Nuclear Norm Minimization,” in Proc. of 2016 IEEE International Symposium on Circuits and Systems (ISCAS 2016), pp. 2639-2642, Montreal, Canada, May 2016. 

