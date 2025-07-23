# pareto.py
"""Multi-objective optimization using Pareto frontier for docking trade-offs."""

def pareto_optimal_docking(cavity_grid, molecule_grid, objectives):
    pareto_frontier = []
    all_solutions = []

    for position in cavity_grid:
        scores = {name: func(cavity_grid, molecule_grid, position) for name, func in objectives.items()}
        all_solutions.append((position, scores))

        is_dominated = False
        for _, other_scores in all_solutions[:-1]:
            if all(other_scores[k] >= scores[k] for k in scores) and any(other_scores[k] > scores[k] for k in scores):
                is_dominated = True
                break

        if not is_dominated:
            pareto_frontier.append((position, scores))

    return pareto_frontier