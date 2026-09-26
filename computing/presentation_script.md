# HorizonRoute presentation

This video uses a synthetic voice. The project used AI assistance.

## 1. HorizonRoute

Horizon Route is a general qubit placement and routing solver for the Q SITE computing challenge. It reduces the public benchmark score from 283 point 5 to 68 point 5. Lower is better. Every output passes the unchanged official checker. The project combines delayed placement with a search over possible routes. This video uses a synthetic voice. OpenAI Codex assisted with implementation, analysis, tests, and writing.

## 2. A route must preserve the program

The input contains logical qubits and a hardware graph. A two qubit operation can run only between connected physical sites. A SWAP moves two logical states and changes where later operations must run. The output must preserve every input operation in its original order, including single qubit operations. The official score adds the SWAP count to one half of the circuit depth. This means that a route must control both movement and time.

## 3. Placement remains open during search

The solver first searches for an embedding of the interaction graph in the hardware graph. A successful embedding needs no SWAPs. Otherwise, the search delays each logical assignment until its first operation. This does not create a state during execution. It assigns a token that was present from the start. After routing, the solver reverses all SWAPs. This recovers a single valid initial placement for the complete program, including states first used late.

## 4. The search considers later operations

For each required interaction, the solver explores shortest routes from either endpoint. It keeps a bounded set of candidates instead of choosing the first available path. One ranking uses future interaction distances. A second ranking also estimates the remaining logical dependency depth. The final choice uses the actual measured score. A complete fallback handles feasible disconnected hardware. The solver contains no benchmark names, saved answers, or lookup table of routes.

## 5. Public benchmark score

Across the six official programs, Horizon Route scores 68 point 5. The supplied baseline scores 283 point 5. This is a reduction of 75 point 84 percent. A stronger comparison method uses static placement search and scores 129 point 5. Horizon Route is 47 point 10 percent lower than that method. The chain and repeated layer cases need zero SWAPs. Each reaches its logical depth lower bound. All routes and measured scores are saved.

## 6. A two-SWAP star route

This figure shows the star program. Logical qubit zero must interact with seven targets. The initial placement puts three targets next to the control. After those interactions, the first SWAP moves the control from site thirteen to site fourteen. Later, a second SWAP moves it to site ten. The control reaches the remaining targets in the correct order. The route uses two SWAPs, has depth nine, and scores six point five.

## 7. Checks and search cost

Ten test groups pass. These include sixty four random programs with mixed operation types. They also test full occupancy, empty programs, unusual labels, disconnected graphs, and invalid inputs. Four additional seed cases score sixty, compared with one hundred sixty one for the official baseline. The default public run took about one hundred one seconds. A smaller search took two point seven nine seconds and scored seventy seven point five. Users can choose this speed and quality tradeoff.

## 8. A complete and repeatable result

The repository contains the required solve function, the unchanged scorer, all measured routes, and an interactive route viewer. The report gives the method, comparisons, sources, and limits. A clean environment can run the tests and reproduce the benchmark table. The bounded search does not prove an optimum for every input. It also does not search every longer detour. Its contribution is a tested general method that sharply reduces the challenge score while preserving the full program.