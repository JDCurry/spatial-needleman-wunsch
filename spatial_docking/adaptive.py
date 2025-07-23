# adaptive.py
"""Scoring function tuning based on known binders and experimental activities."""

from scipy.optimize import minimize
import numpy as np

class AdaptiveScoringFunction:
    def __init__(self, base_weights):
        self.weights = base_weights

    def update_weights(self, cavity_grid, known_binders, known_activities, scoring_function):
        def loss(w):
            self.weights = dict(zip(self.weights.keys(), w))
            predicted = [scoring_function(cavity_grid, b, self.weights) for b in known_binders]
            corr = np.corrcoef(predicted, known_activities)[0, 1]
            return -corr  # maximize correlation

        initial = list(self.weights.values())
        result = minimize(loss, initial, method='Nelder-Mead')
        self.weights = dict(zip(self.weights.keys(), result.x))
        return self.weights