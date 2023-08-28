Script used to perform differential expression on data mini-bulk expression data.

The test performed is a variance adjusted t-test between compound treated conditions and negative control conditions. Control samples can be randomized to reduce systemic biases associated with using the same control group for each comparison.

The variance adjustment controls for low replicate numbers in the treatment conditions by first fitting a curve between the mean expression and coefficient of variation for each treatment. Condtions with lower variance than predicted by the regression will be adjusted to the predicted value

### Example

`python3 variance_adjusted_ttest_modified.py -i new_rpkms.txt -o modified_script -r New -v linear -m meta.txt -l 1hr -c 1hr --subsample_g1 10 --test_dict_path break_dict_for_DE_new.json --replicate_column compound_name`
