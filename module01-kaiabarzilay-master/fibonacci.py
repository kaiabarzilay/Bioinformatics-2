#!/usr/bin/env python
# fibonacci.py

def population(n, k): # n = days, k = rate
    """
    Find the population on a certain day (n), based on the rate of reproduction (k)
    """
    population_value = [0] * (n + 1)
    # For example population at day 1 is 1
    population_value[0] = 0
    population_value[1] = 1
    for current_value in range(2, n+1):
        # Equation to find population on a given day
        population_value[current_value] = k * (population_value[current_value-2]) + population_value[current_value-1]
    return population_value[n]

# Example output for population(4,6) is 13
# final_population = population(4, 6)
# print(final_population)
