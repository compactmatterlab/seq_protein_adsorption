# seq_protein_adsorption

How to use - 

All files are MATLAB files, compatible with MATLAB version 2016B and later. 

To compare model predictions with experimental data, you can run the following command in MATLAB command window:

loglog(time_series_data,survival fraction_data,'ko',time_series_data,myadsorption_gillespe(time_series_data,survival fraction_data,tau,p,N),'b-')

or 

loglog(time_series_data,survival fraction_data,'ko',time_series_data,myadsorption_gillespe_rand(time_series_data,survival fraction_data,tau,p,N),'b-')

myadsorption_gillespe assumes same p and 1-p values for all steps, myadsorption_gillespe_rand randomly generates unique p and 1-p values for each step, normally distributed around mean p and 0.025 standard deviation. 

Alternative, you can use the myadsorption_gillespe_code and myadsorption_gillespe_rand_code files instead of the above mentioned function files. 

For finding the best fit parameters to the model, use the ParameterOptimizer_loglog.m file. To be able to generate a good fit within reasonable time, you will have to start with a good initial guess of the paramaters. Parameter fitting is obtained by minimizing the sum of squares between log(experimental) and log(predicted) survival fractions and performing a Metropolis Algorithm based restricted Monte Carlo search of the parameter space. 

For any question or further details, please email pkatira@sdsu.edu
