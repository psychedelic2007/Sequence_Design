from deap import base, creator, tools, algorithms
import random
import numpy as np
import pandas as pd

df = pd.read_csv("new_results_13kseqs/variants/final_results/Scores/merged_scores.csv")
creator.create("FitnessMulti", base.Fitness, weights=(1.0, -1.0))
creator.create("Individual", list, fitness=creator.FitnessMulti)

toolbox = base.Toolbox()

#[wm, we, wb, wt], each in [0.0, 1.0]
toolbox.register("attr_float", random.uniform, 0.0, 1.0)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, 4)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

def evaluate(individual):
    wm, we, wb, wt = individual
    total = wm + we + wb + wt
    if total == 0:
        total = 1e-10
    wm, we, wb, wt = wm / total, we / total, wb / total, wt / total

    modify_flags = []
    for _, row in df.iterrows():
        mut, esc, bc, tc = (
            row["mutation_frequency"],
            row["escape"],
            row["bcell_score"],
            row["tcell_score"],
        )
        score = (wm * mut + we * esc) - (wb * bc + wt * tc)
        modify_flags.append(1 if score > 0 else 0)

    total_escape = sum(
        (row["escape"] + row["mutation_frequency"])
        for j, (_, row) in enumerate(df.iterrows())
        if modify_flags[j] == 1
    )

    epitope_loss = sum(
        1
        for j, (_, row) in enumerate(df.iterrows())
        if modify_flags[j] == 1
        and (row["bcell_score"] >= 0.9 or row["tcell_score"] >= 0.9)
    )

    return total_escape, epitope_loss


toolbox.register("evaluate", evaluate)
toolbox.register("mate", tools.cxBlend, alpha=0.5)
toolbox.register("mutate", tools.mutPolynomialBounded, eta=20.0, low=0.0, up=1.0, indpb=0.3)
toolbox.register("select", tools.selNSGA2)

pop = toolbox.population(n=100)
hof = tools.ParetoFront()
stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg", np.mean, axis=0)
stats.register("min", np.min, axis=0)
stats.register("max", np.max, axis=0)

algorithms.eaMuPlusLambda(
    pop,
    toolbox,
    mu=100,
    lambda_=200,
    cxpb=0.5,
    mutpb=0.3,
    ngen=40,
    stats=stats,
    halloffame=hof,
    verbose=True,
)

pareto_results = pd.DataFrame(
    [
        {
            "mutation_weight": ind[0] / sum(ind) if sum(ind) > 0 else 0,
            "escape_weight": ind[1] / sum(ind) if sum(ind) > 0 else 0,
            "bcell_weight": ind[2] / sum(ind) if sum(ind) > 0 else 0,
            "tcell_weight": ind[3] / sum(ind) if sum(ind) > 0 else 0,
            "escape_removed": ind.fitness.values[0],
            "epitope_loss": ind.fitness.values[1],
        }
        for ind in hof
    ]
)

pareto_results.to_csv("pareto_frontier_results.csv", index=False)
print("Optimization complete. Pareto front saved to pareto_frontier_results.csv")
