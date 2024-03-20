import numpy as np
def mean_bias_deviation(y_true, y_pred):
    """
    Compute the Mean Bias Deviation (MBD)
    Positive values indicate underestimation on average, 
    while negative values indicate overestimation.
    """
    return np.mean((y_true - y_pred) / y_true)