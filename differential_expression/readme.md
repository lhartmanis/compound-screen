Script used to perform differential expression on data mini-bulk expression data.

The test performed is a variance adjusted t-test between compound treated conditions and negative control conditions. Control samples can be randomized to reduce systemic biases associated with using the same control group for each comparison.

The variance adjustment controls for low replicate numbers in the treatment conditions by first fitting a curve between the mean expression and coefficient of variation for each treatment. Condtions with lower variance than predicted by the regression will be adjusted to the predicted value


### Insert example of how to run this analysis
