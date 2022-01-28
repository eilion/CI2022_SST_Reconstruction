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