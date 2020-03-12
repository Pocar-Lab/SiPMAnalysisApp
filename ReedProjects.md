# Projects Questions and Progress

## Analysis Procedures and the VUV4 Dataset:

- In calculating Chi^2 we want to know the variance of each bin. Intuitively a bin follows a binomial distribution from run to run but most procedures on the internet suggest using E_i, which would work if the bin followed a Poissonian distribution.

- Why does the Scipy curve_fit covariance matrix yield identical values for the error on the estimators of the mean and sigma? My progress so far indicates that they should vary by a factor of sqrt(2)

- Michelle wants me to add a feature by which you can plot multiple lines in the plotting applet

- Analysis fits the entire dataset but only uses a limited section to calculate Chi^2. I should only fit the limited section.

- Andrea wants to see that the error on the estimator of the mean is a factor of 1/sqrt(N) different from the estimator of the standard deviation

- Eventually the table should include a collumn either giving the number of datapoints eventually used or the final value for the error on the midpoint.

- Continue to work on the write-up

### Why the denominator of the $\chi^2$ statistic is equal to $\alpha_i^2$ in Measurements and their Uncertainties textbook but is equal to $E_i^2$ in every other source?

__Other Sources:__
- https://mathworld.wolfram.com/Chi-SquaredTest.html
- https://en.wikipedia.org/wiki/Chi-squared_distribution

* It seems that minimization of $\chi^2$ is a particular case of a weighted least squares minimization: https://en.wikipedia.org/wiki/Weighted_least_squares

- "The weights should, ideally, be equal to the reciprocal of the variance of the measurement."

* Eventually it is the case that minimizing least squares is equivelant to the maximizing the Likelihood function: https://www.tandfonline.com/doi/pdf/10.1080/01621459.1976.10481508?needAccess=true

* Out of curiousity I calculated the Likelihood function for a fit to a histogram. Which is as such:
![L = \prod_{i} {N}\choose{O_i} \left(\frac{E_i}{N}\right)^N \left(1-\frac{E_i}{N}\right)^N-O_i](https://render.githubusercontent.com/render/math?math=L%20%3D%20%5Cprod_%7Bi%7D%20%7BN%7D%5Cchoose%7BO_i%7D%20%5Cleft(%5Cfrac%7BE_i%7D%7BN%7D%5Cright)%5EN%20%5Cleft(1-%5Cfrac%7BE_i%7D%7BN%7D%5Cright)%5EN-O_i)

__Thought Process of the derivation:__
1) We note that the height of a bin is a binomial random variable
2) The probability to obtain the observed height O_i is *binom(O_i, p = E_i/N, N)*
3) The probability to obtain the *exact* arrangement of bins with the observed heights is the product of the probability to obtain each bins observed height over all the bins


### The Covariance matrix pcov returned by sp.curve_fit():

- The source code for curve_fit seems to return directly the value cov_x obtained from the function leastsq: 
res = leastsq(func, p0, Dfun=jac, full_output=1, kwargs)
popt, pcov, infodict, errmsg, ier = res

- The documentation for leastsq function says:
cov_x : The inverse of the Hessian. fjac and ipvt are used to construct an estimate of the Hessian. A value of None indicates a singular matrix, which means the curvature in parameters x is numerically flat. To obtain the covariance matrix of the parameters x, cov_x must be multiplied by the variance of the residuals – see curve_fit."
- It doesn't appear that this suggested step occurs. Though it's unclear to me what the "variance of the residuals" is exactly. EDIT: It does happen, this is the code:
elif not absolute_sigma:
        if ysize > p0.size:
            s_sq = cost / (ysize - p0.size)
            pcov = pcov * s_sq
- Another line says cost is the sum of squares, and ysize - p0.size is simply N - #params, aka s_sq is the minimized chi squared, as prescribed.
- This process only happens if the errors are not provided (aka assumed homogeneous). I am now following this article to determine the effect of Sigma. https://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-using-the-optimize-leastsq-method-i

"optimize.leastsq and optimize.curvefit provide us a way to estimate errors in fitted parameters, but we cannot just use these methods without questioning them a little bit. The bootstrap is a statistical method which uses brute force, and in my opinion, it has a tendency of working better in situations that may be harder to interpret."

- Ultimately this is all just a wrapper around MINPACK algorithms which are written in FORTRAN, there is no way I am going to understand the code, even the approximation to the Hessian matrix is going to be too challenging.
- If I am to trust curve_fit documentation, this is what it says: The estimated covariance of popt. The diagonals provide the variance of the parameter estimate. To compute one standard deviation errors on the parameters use perr = np.sqrt(np.diag(pcov)).
