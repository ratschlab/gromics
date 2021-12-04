import numpy as np

def total_count(expression):
    """
    Compute total count library size.

    :param expression: Expression count matrix with genes as rows and samples as columns
    :return: integer total count library size per sample
    """
    return expression.sum(axis=0) 

def upper_quartile(expression):
    """
    Compute total count library size.

    :param expression: Expression count matrix with genes as rows and samples as columns
    :return: float upper quartile of expression per sample
    """
    return np.array([np.percentile(x, 75) for x in expression.T])
