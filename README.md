
# ProJect

## Introduction

ProJect is a novel mixed-model imputation method, as a powerful and versatile solution for handling MVs in diverse data modalities, including genomics and mass-spectrometry-based proteomics. ProJect demonstrates consistent and superior performance as evidenced by its success on various benchmark datasets.

## Usage

### Function: `initial(data2, subdata2, flag)`

**Parameters:**

- `data2`: The column of data that you want to impute. We recommend performing a log transformation on this column before imputation.
- `subdata2`: The rest of the data excluding the column you want to impute. We also recommend log transformation for this dataset before imputation.
- `flag`: 
  - A value of `-1` means **pure MNAR (Missing Not at Random) imputation**.
  - Any other value will trigger **mixed imputation** (a combination of MNAR and MAR).

**Note**: If you wish to impute the entire dataset column by column, you will need to use a for loop and call the function iteratively for each column.

### Example:

Here’s an example of how you can use this function to impute the first column of a dataset named `RC`:

```r
# Prepare the data
RCN1 = log(RC[, 1])  # Take the first column and apply log transformation
restRC1 = log(RC[, -1])  # Apply log transformation to the rest of the dataset

# Perform the imputation
impRC = initial(RCN1, restRC1, -1)  # Start the imputation process with MNAR flag (-1)

# Extract the results
impdata = impRC$imp.data  # This is the normalized output data
meand = impRC$meand  # This is the mean value of the output
sdd = impRC$sdd  # This is the standard deviation of the output

# Reverse the log transformation to get the final imputed column
imptmp = exp(impdata * sdd + meand)  # Final imputed first column
```

For further details, refer to the `Project_func.R` file for additional functions and examples. For the full script, containing benchmark and comparisons, please refer to the `project.Rmd` file.

