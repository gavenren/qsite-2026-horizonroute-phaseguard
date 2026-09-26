# HorizonRoute demonstration

HorizonRoute solves the placement and routing problem together. Its public score is sixty-eight point five. The official baseline scores two hundred eighty-three point five. This is a seventy-six percent reduction, with no change to the scorer.

The task has a simple cost rule. Each SWAP costs one point. Each circuit layer costs half a point. A valid result must execute every original operation in the same order. It must also place every two-qubit operation on a hardware edge.

A poor initial placement makes this task much harder. The baseline fixes each logical qubit at the start. It then moves one endpoint along a shortest path. It does not consider how that movement affects later gates.

Our main idea is delayed placement. During search, an unused physical position holds an unassigned token. When a logical qubit first appears, the solver chooses a token for it. At the end, the solver reverses every SWAP. This finds where that token started. The result is a normal initial placement for the complete circuit.

The solver first checks whether the interaction graph fits the hardware without any SWAPs. This solves the chain and the repeated-layer case directly. For harder inputs, it keeps several candidate routes after each operation.

Each candidate records the qubit positions, the inserted SWAPs, and the last-used layer on every physical qubit. The search considers movement from both endpoints. It estimates the cost of later interactions. A second pass also estimates the remaining logical critical path. The solver returns the candidate with the best actual score.

Now look at the star example in the viewer. Logical qubit zero must meet seven different partners. The hardware degree is at most three. The route starts with nearby partners. Two SWAPs move the center to reach the remaining partners. The final result has two SWAPs and nine layers. Its score is six point five.

That star score meets a lower bound. One center SWAP can expose at most two new partners. At least two are necessary. With exactly two, both must move the center, so the seven center gates and two SWAPs need nine layers. Three SWAPs with seven layers have the same minimum score. The chain and repeated-layer results also meet their lower bounds.

The harder cases show why joint search matters. The dense case falls from one hundred twenty-two points to thirty-five point five. A separate static-placement search scores seventy-four on that case. On all six public cases, static placement scores one hundred twenty-nine point five. Our result is forty-seven percent lower.

We checked more than the public score. Ten test groups pass, including sixty-four random mixed-operation cases. Four additional random circuits score sixty points in total, compared with one hundred sixty-one for the official baseline. The original scorer checks every result.

The full public run took about one hundred one seconds. A smaller search took under three seconds and scored seventy-seven point five. The package includes both settings, complete routes, test results, and a local viewer.

This is a bounded search. We do not claim a global optimum for every circuit. The implementation uses no stored benchmark answers. OpenAI Codex assisted with the work. The included scripts reproduce the measured evidence.
