# HorizonRoute: Joint Placement and Qubit Routing

QSITE 2026 · Computational Track

## Result

HorizonRoute scores **68.5** on the six public cases. The official baseline scores **283.5**. The score decreases by **75.84%**. Every output passes the unchanged official scorer. Lower scores are better.

The solver also beats a static-placement search. That comparison uses simulated annealing, then greedy shortest-path routing. It scores 129.5. HorizonRoute reduces this score by 47.10%.

| Public case | Official baseline | Static placement | HorizonRoute | SWAPs | Depth |
|---|---:|---:|---:|---:|---:|
| ghz_star | 14.0 | 11.0 | 6.5 | 2 | 9 |
| chain_trotter | 15.0 | 4.5 | 4.5 | 0 | 9 |
| ladder_trotter | 35.5 | 10.0 | 6.5 | 3 | 7 |
| qaoa_random | 39.0 | 21.5 | 12.5 | 6 | 13 |
| dense_random | 122.0 | 74.0 | 35.5 | 24 | 23 |
| vqe_layers | 58.0 | 8.5 | 3.0 | 0 | 6 |
| **Total** | **283.5** | **129.5** | **68.5** | **35** | **67** |

## Main idea

A fixed initial placement can force unnecessary movement. HorizonRoute delays each placement decision until the logical qubit first appears. It searches placement and routing together.

This method does not create a state during execution. Each unused physical position holds an unassigned token. A first-use decision assigns a logical state to that token. At the end, the solver reverses every SWAP. This gives the token's initial physical position. Distinct tokens give distinct initial positions.

The solver first searches for a placement with no SWAPs. If it finds one, it returns that placement directly. Otherwise, it keeps a beam of candidate routes after each input operation.

For the next gate, it considers shortest routes from either endpoint. Each candidate keeps its placement, SWAP count, operation list, and last-used layer on every physical qubit. Its actual cost is the official objective: SWAP count plus half the depth.

Two search passes use different future estimates. One estimates distances for the next 16 gates. The other also estimates the remaining logical critical path. The passes keep up to 900 and 500 states. The solver selects their best actual result. It uses no benchmark names or stored answers.

---PAGE---

## Correctness and validation

The returned list preserves every input operation, its operands, and its position in the sequence. The solver only inserts SWAPs. Every two-qubit operation uses a hardware edge. Each SWAP updates the logical-to-physical map.

The official scorer then translates the entire result back to logical qubits. Its output must equal the original input exactly. The scorer also assigns parallel layers with its original scheduling function.

Ten test groups pass. They include 64 random cases, mixed single-qubit operations, non-contiguous labels, repeated pairs, and full occupancy. Tests also cover empty input, disconnected hardware, impossible capacity, and deterministic output.

Four additional cases use seeds 4101–4104. Each has eight logical qubits, twenty random two-qubit gates, and one single-qubit gate. These cases were not used to choose the final settings. HorizonRoute scores 60.0, versus 161.0 for the official baseline and 83.5 for static placement. These are small local tests, not the organizer's hidden tests.

## Three results meet lower bounds

The chain needs nine ordered layers. HorizonRoute uses no SWAPs and scores 4.5. The repeated-layer case needs six layers. Its score of 3.0 is also optimal.

The star has seven distinct partners. Hardware degree is at most three. A center SWAP can expose at most two new partners. Another SWAP can expose at most one. Thus, at least two SWAPs are necessary. With exactly two, both must move the center, which requires at least nine layers. With three or more, the seven center gates still require seven layers. Both bounds give a score of at least 6.5. HorizonRoute attains it.

## Cost and limits

The full public run took 100.82 seconds on the recorded local system. A beam width of 40 took 2.79 seconds and scored 77.5. This gives a useful speed option.

The search is bounded. It keeps at most 96 shortest routing sequences for a physical pair. It can miss useful longer paths. We claim no optimum for the ladder or random cases. We also claim no quantum hardware result or win probability.

## Sources and disclosure

The unchanged challenge files are from Quantum Coalition QSITE 2026, commit `a45b04211d4da27dca3c79dbe7009ce58e85dc48`. The package includes their MIT notice. Source: https://github.com/benmcdonough20/QSITE-2026-QuantumCoalition

SABRE motivates look-ahead routing; this solver uses a different search method. G. Li, Y. Ding, and Y. Xie, ASPLOS 2019: https://arxiv.org/abs/1809.02573

OpenAI Codex assisted with code, tests, analysis, and writing. Executed scripts produced all reported scores.
