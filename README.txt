To run the software, you should follows the below steps:
1. Make a folder 'Inputs/INPUT'.
2. Fill in the folder by following the instruction in 'Inputs/README'.
3. Define 'Inputs/INPUT/setting.txt' as follows:
3.1. 'number_of_samples' is the number of posterior sea surface temperature samples to draw. Default is 10000.
3.2. 'max_iters' is the number of iterations to estimate kernel hyperparameters. Default is 10000.
3.3. 'mean_T' and 'stdv_T' define the shift and scale hyperparameters to standardize sea surface temperatures. Defaults are both 'NaN', which means that these two values are automatically determined based on the input records.
3.4. 'query_age_start' and 'query_age_end' define the start and end query ages of the SST inference. Defaults are both 'NaN', which means that these two values are automatically determined based on the input records.
3.5. 'number_of_query_ages' defines how many query ages will be defined between 'query_age_start' and 'query_age_end'. Default is 101.
4. Define 'Inputs/INPUT/MEAN.txt' by the prior SST mean. The first column indicates ages and second one defines the associated SSTs (in degree in celsius).
5. Run the software by following the instruction in 'main.m'.

Once the software terminates the steps, you will get results in 'Outputs/INPUT_KERNEL'. 'result.mat' consists of four structures.
1. 'data' summarizes the results.
1.1. 'data.name' indicates the name of an input record.
1.2. 'data.T' stores the ages.
1.3. 'data.Y' stores the UK'37 records.
1.4. 'data.T0' stores the query ages.
1.5. 'data.MEAN0' stores the prior SST mean values at 'data.T0' that is defined by 'Inputs/INPUT/MEAN.txt'.
1.6. 'data.MEAN' stores the prior SST mean values at 'data.T' that is defined by 'Inputs/INPUT/MEAN.txt'.
1.7. 'data.SCALE' and 'data.BIAS' store the estimated core-specific scale and shift parameters.
1.8. 'data.LB_INDIV', 'data.MD_INDIV' and 'data.UB_INDIV' store the 2.5%, 50% and 97.5% quantiles of the posterior SST samples at 'data.T' by the point-wise SST reconstruction.
1.9. 'data.LB', 'data.MD' and 'data.UB' store the 2.5%, 50% and 97.5% quantiles of the posterior SST samples at 'data.T' by the GPST SST reconstruction.
1.10. 'data.LB0', 'data.MD0' and 'data.UB0' store the 2.5%, 50% and 97.5% quantiles of the posterior SST samples at 'data.T0' by the GPST SST reconstruction.
2. 'param' stores the parameters of Gaussian process transition model.
2.1. 'param.GAMMA' and 'param.LAMBDA' store the estimated kernel hyperparameters.
2.2. 'param.LOGLIK' stores the values of objective function over iterations.
3. 'setting' stores the settings defined in 'Defaults/setting.txt' and 'Inputs/INPUT/setting.txt'.
4. 'Samples' stores the posterior SST samples.