# HorizonRoute

A general qubit router for the QSITE 2026 Computational Track.

The solver chooses an initial placement. It inserts legal SWAP operations and preserves the exact input operation order. It uses no saved benchmark answers.

**Measured public score: 68.5**, compared with 283.5 for the supplied baseline and 129.5 for static placement search. Lower is better. All 10 test groups pass.

Read the [two-page report](report.pdf), view the [presentation slides](presentation-reviewed.pptx), or download the [complete source package](computing-source.zip). The package preserves the required data and starter-code folders.

## Start

Use Python 3.14 for the same environment as the challenge. The core solver uses the Python standard library and a NetworkX graph.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python test_solver.py
python benchmark.py
python make_figures.py
```

Open `demo.html` in a browser for the route viewer. No server is necessary.

Use the required interface:

```python
from solver import solve
placement, routed_program = solve(program, hardware_graph)
```

For a quick result, use `solve(program, hardware_graph, beam_width=40, portfolio=False)`. The default uses two search methods. Their beam widths are 900 and 500. More search can improve a result, but it does not guarantee an improvement.

## Method

1. Search for an initial placement that needs no SWAPs.
2. If this fails, delay each logical assignment until its first operation.
3. Keep several candidate routes after each input operation.
4. Search shortest routes from both endpoints of the next required gate.
5. Rank candidates with the measured prefix cost and estimated future costs.
6. Compare a distance estimate with an estimate that also uses the remaining logical critical path.
7. Select the route with the lowest actual score.
8. Reverse its SWAP sequence to recover the initial placement.

The delayed assignment is a search method. It does not prepare a qubit during execution. Each logical state occupies a distinct token from the start. The final reverse pass finds that token's initial physical position.

A fallback packs logical components into hardware components. It gives a valid route for feasible disconnected hardware. The beam search is used on connected hardware.

## Files

- `solver.py`: required `solve()` function.
- `benchmark.py`: public cases, a random-placement baseline, an ablation, and unseen cases.
- `test_solver.py`: correctness and error tests with the unchanged official scorer.
- `results/benchmark_results.json`: all placements, routes, layers, scores, and run times.
- `results/scores.csv`: compact public results.
- `report.md`: the challenge report.
- `demo_script.md`: narration for the demonstration.
- `demo.html`: local route viewer.
- `figures/`: measured score and route figures.
- `vendor/starter_kit/`: unchanged official starter files.

## Scope and limits

The search has a fixed width. It does not prove an optimum for every input. It considers up to 96 shortest routing sequences for each physical pair. It does not search every longer detour. The two zero-SWAP results meet their logical depth lower bounds. The star result also meets a simple lower bound; see the report.

The router supports `1Q` and `2Q` operations. The challenge represents these as abstract operations. The router does not cancel gates or change their physical definitions. Such changes would violate the scorer's exact operation-order check.

The public cases were used during development. The unseen cases use separate seeds. They are not the organizer's hidden test set. They do not establish a win probability.

## Sources and disclosure

The challenge files are from commit `a45b04211d4da27dca3c79dbe7009ce58e85dc48` of [Quantum Coalition QSITE 2026](https://github.com/benmcdonough20/QSITE-2026-QuantumCoalition). Their MIT notice is in `vendor/LICENSE`. All six Python starter modules are unchanged.

[SABRE, Li, Ding, and Xie, ASPLOS 2019](https://arxiv.org/abs/1809.02573) is a source for look-ahead routing and the placement problem. HorizonRoute uses a different search method with strict input order. It is not an implementation of SABRE.

OpenAI Codex assisted with code, tests, analysis, and writing. Executed test scripts produced the reported results. No result is a quantum hardware measurement.
