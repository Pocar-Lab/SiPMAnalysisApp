# Projects Questions and Progress

## Analysis Procedures and the VUV4 Dataset:

- In calculating Chi^2 we want to know the variance of each bin. Intuitively a bin follows a binomial distribution from run to run but most procedures on the internet suggest using E_i, which would work if the bin followed a Poissonian distribution.

- Why does the Scipy curve_fit covariance matrix yield identical values for the error on the estimators of the mean and sigma? My progress so far indicates that they should vary by a factor of sqrt(2)

- Michelle wants me to add a feature by which you can plot multiple lines in the plotting applet

- Analysis fits the entire dataset but only uses a limited section to calculate Chi^2. I should only fit the limited section.

- Andrea wants to see that the error on the estimator of the mean is a factor of 1/sqrt(N) different from the estimator of the standard deviation

- Eventually the table should include a collumn either giving the number of datapoints eventually used or the final value for the error on the midpoint.

- Continue to work on the write-up
