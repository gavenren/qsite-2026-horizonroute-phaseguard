"""HorizonRoute: general qubit placement and routing with bounded beam search.

Only the Python standard library and the NetworkX graph interface are needed.
All input operations remain in the same order. The solver does not use benchmark
names, stored answers, scorer internals, or a database of circuits.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from functools import lru_cache
from heapq import nsmallest
from math import inf


@dataclass(slots=True)
class _State:
    positions: tuple[int, ...]
    last: tuple[int, ...]
    swaps: int
    route: tuple[tuple, ...]

    @property
    def cost(self):
        return self.swaps + 0.5 * max(self.last, default=0)


def _embedding(program, adjacency, m, node_budget=50000):
    """Find a zero-SWAP graph embedding before the more costly route search."""
    neighbors = [set() for _ in range(m)]
    for op in program:
        if op[0] == '2Q':
            neighbors[op[1]].add(op[2])
            neighbors[op[2]].add(op[1])
    if max(map(len, neighbors), default=0) > max(map(len, adjacency), default=0):
        return None
    placement = [-1] * m
    used = set()
    visits = 0

    def visit():
        nonlocal visits
        visits += 1
        if visits > node_budget:
            return None
        remaining = [q for q in range(m) if placement[q] < 0]
        if not remaining:
            return tuple(placement)
        logical = max(remaining, key=lambda q: (sum(placement[k] >= 0 for k in neighbors[q]), len(neighbors[q]), -q))
        linked = [placement[k] for k in neighbors[logical] if placement[k] >= 0]
        candidates = set(range(len(adjacency))) - used
        for p in linked:
            candidates.intersection_update(adjacency[p])
        candidates = sorted(candidates, key=lambda p: (abs(len(adjacency[p]) - len(neighbors[logical])), p))
        for physical in candidates:
            if len(adjacency[physical]) < len(neighbors[logical]):
                continue
            placement[logical] = physical
            used.add(physical)
            result = visit()
            if result is not None:
                return result
            used.remove(physical)
            placement[logical] = -1
        return None
    return visit()


def _search(program, adjacency, distances, m, width, lookahead_weight, initial_positions=None, critical=True):
    n = len(adjacency)

    @lru_cache(maxsize=None)
    def moves(left, right):
        """All shortest ways to make the next pair adjacent, from either end."""
        if distances[left][right] == 1:
            return ((),)
        if distances[left][right] == inf:
            return ()
        result = []
        for next_left in adjacency[left]:
            if distances[next_left][right] == distances[left][right] - 1:
                for tail in moves(next_left, right):
                    result.append(((left, next_left),) + tail)
        for next_right in adjacency[right]:
            if distances[left][next_right] == distances[left][right] - 1:
                for tail in moves(left, next_right):
                    result.append(((right, next_right),) + tail)
        # Equivalent interleavings can produce the same mapping. A fixed cap
        # controls memory on graphs with very many equal shortest paths.
        return tuple(dict.fromkeys(result))[:96]

    future = []
    for index in range(len(program)):
        pairs = [(op[1], op[2]) for op in program[index + 1:] if op[0] == '2Q'][:16]
        future.append([(a, b, 0.90 ** j) for j, (a, b) in enumerate(pairs)])

    def rank(state, index):
        pos = state.positions
        used = set(pos)
        free = tuple(p for p in range(n) if p not in used)
        forecast = 0.0
        unknown = {}
        for a, b, weight in future[index]:
            pa, pb = pos[a], pos[b]
            if pa >= 0 and pb >= 0:
                forecast += weight * max(0, distances[pa][pb] - 1)
            elif pa >= 0 or pb >= 0:
                fixed = max(pa, pb)
                if critical:
                    latent = b if pa >= 0 else a
                    unknown.setdefault(latent, []).append((fixed, weight))
                else:
                    forecast += weight * min((max(0, distances[fixed][p] - 1) for p in free), default=0)
        if not critical:
            return state.cost + lookahead_weight * forecast
        for terms in unknown.values():
            forecast += min((sum(weight * max(0, distances[fixed][p] - 1)
                                 for fixed, weight in terms) for p in free), default=0)
        # Estimate the remaining depth from logical dependencies. This avoids
        # choosing an early shallow prefix that hides a long later path.
        logical_last = [state.last[p] if p >= 0 else 0 for p in pos]
        for gate in program[index + 1:]:
            if gate[0] == '2Q':
                a, b = gate[1:]
                layer = 1 + max(logical_last[a], logical_last[b])
                logical_last[a] = logical_last[b] = layer
        projected_depth = max(max(state.last), max(logical_last))
        return state.swaps + 0.5 * projected_depth + lookahead_weight * forecast

    states = [_State(tuple(initial_positions) if initial_positions is not None else (-1,) * m, (0,) * n, 0, ())]
    for index, op in enumerate(program):
        next_states = {}
        for state in states:
            pos = state.positions
            free = [p for p in range(n) if p not in pos]
            a = op[1]
            assignments = []
            if op[0] == '1Q':
                for pa in ([pos[a]] if pos[a] >= 0 else free):
                    newpos = list(pos)
                    newpos[a] = pa
                    new = _State(tuple(newpos), state.last, state.swaps, state.route + (('1Q', pa),))
                    key = (new.positions, new.last)
                    old = next_states.get(key)
                    if old is None or new.swaps < old.swaps:
                        next_states[key] = new
                continue
            b = op[2]
            free_has_edge = any(q in free for p in free for q in adjacency[p])
            lefts = [pos[a]] if pos[a] >= 0 else free
            rights = [pos[b]] if pos[b] >= 0 else free
            for pa in lefts:
                for pb in rights:
                    if pa == pb:
                        continue
                    # With two new logical states, an adjacent starting pair is
                    # sufficient for this bounded search. Later swaps are legal.
                    if pos[a] < 0 and pos[b] < 0 and free_has_edge and distances[pa][pb] != 1:
                        continue
                    assignments.append((pa, pb))
            for pa, pb in assignments:
                assigned = list(pos)
                assigned[a], assigned[b] = pa, pb
                for steps in moves(pa, pb):
                    newpos = assigned.copy()
                    last = list(state.last)
                    extra = []
                    for left, right in steps:
                        for q, p in enumerate(newpos):
                            if p == left:
                                newpos[q] = right
                            elif p == right:
                                newpos[q] = left
                        layer = 1 + max(last[left], last[right])
                        last[left] = last[right] = layer
                        extra.append(('SWAP', left, right))
                    left, right = newpos[a], newpos[b]
                    layer = 1 + max(last[left], last[right])
                    last[left] = last[right] = layer
                    extra.append(('2Q', left, right))
                    new = _State(tuple(newpos), tuple(last), state.swaps + len(steps), state.route + tuple(extra))
                    key = (new.positions, new.last)
                    old = next_states.get(key)
                    if old is None or new.swaps < old.swaps:
                        next_states[key] = new
        if not next_states:
            raise ValueError('The search found no route. Check graph connectivity and qubit capacity.')
        states = nsmallest(width, next_states.values(), key=lambda s: (rank(s, index), s.cost, s.positions, s.last))
    return min(states, key=lambda s: (s.cost, s.swaps, s.positions))



def _fallback(program, adjacency, distances, m):
    """Make a valid route even if the hardware has disconnected components."""
    def components(edges):
        unseen = set(range(len(edges)))
        groups = []
        while unseen:
            seed = min(unseen)
            unseen.remove(seed)
            group, queue = [seed], [seed]
            while queue:
                p = queue.pop()
                for q in edges[p]:
                    if q in unseen:
                        unseen.remove(q)
                        group.append(q)
                        queue.append(q)
            groups.append(sorted(group))
        return sorted(groups, key=lambda group: (-len(group), group))

    interaction = [set() for _ in range(m)]
    for op in program:
        if op[0] == '2Q':
            interaction[op[1]].add(op[2])
            interaction[op[2]].add(op[1])
    logical_groups = components(interaction)
    physical_groups = components(adjacency)
    capacity = [len(group) for group in physical_groups]
    assignment = []

    def pack(index):
        if index == len(logical_groups):
            return True
        size = len(logical_groups[index])
        tried = set()
        for component, room in enumerate(capacity):
            if room >= size and room not in tried:
                tried.add(room)
                capacity[component] -= size
                assignment.append(component)
                if pack(index + 1):
                    return True
                assignment.pop()
                capacity[component] += size
        return False
    if not pack(0):
        raise ValueError('Hardware components cannot hold the connected logical components.')
    pools = [list(group) for group in physical_groups]
    pos = [-1] * m
    for group, component in zip(logical_groups, assignment):
        for q in group:
            pos[q] = pools[component].pop(0)
    last = [0] * len(adjacency)
    route, swaps = [], 0
    for op in program:
        if op[0] == '1Q':
            route.append(('1Q', pos[op[1]]))
            continue
        a, b = op[1:]
        while distances[pos[a]][pos[b]] > 1:
            left = pos[a]
            right = next(p for p in adjacency[left] if distances[p][pos[b]] == distances[left][pos[b]] - 1)
            pos = [right if p == left else left if p == right else p for p in pos]
            layer = 1 + max(last[left], last[right])
            last[left] = last[right] = layer
            route.append(('SWAP', left, right))
            swaps += 1
        left, right = pos[a], pos[b]
        layer = 1 + max(last[left], last[right])
        last[left] = last[right] = layer
        route.append(('2Q', left, right))
    return _State(tuple(pos), tuple(last), swaps, tuple(route))


def solve(program, hardware_graph, *, beam_width=900, lookahead_weight=0.70, portfolio=True):
    """Return (initial_placement, routed_program) in the official interface.

    Optional search parameters do not change the two-argument interface.
    More beam width uses more time and can improve the score. The result is
    deterministic for the same input, graph, and settings.
    """
    if not isinstance(beam_width, int) or beam_width < 1:
        raise ValueError('beam_width must be a positive integer.')
    if hardware_graph.is_directed():
        raise ValueError('Use an undirected hardware graph.')
    program = [tuple(op) for op in program]
    nodes = sorted(hardware_graph.nodes, key=lambda p: (str(type(p)), repr(p)))
    logical = sorted({q for op in program for q in op[1:]})
    if len(logical) > len(nodes):
        raise ValueError('There are more logical qubits than physical qubits.')
    if not program:
        return {}, []
    for op in program:
        if op[0] not in ('1Q', '2Q') or len(op) != (2 if op[0] == '1Q' else 3):
            raise ValueError(f'Unsupported operation: {op!r}')
        if op[0] == '2Q' and op[1] == op[2]:
            raise ValueError('A two-qubit gate must use two distinct logical qubits.')
    physical_index = {p: i for i, p in enumerate(nodes)}
    logical_index = {q: i for i, q in enumerate(logical)}
    internal = [(op[0], *(logical_index[q] for q in op[1:])) for op in program]
    adjacency = tuple(tuple(sorted(physical_index[q] for q in hardware_graph.neighbors(p))) for p in nodes)
    n, m = len(nodes), len(logical)
    distances = []
    for source in range(n):
        row = [inf] * n
        row[source] = 0
        queue = deque([source])
        while queue:
            p = queue.popleft()
            for q in adjacency[p]:
                if row[q] == inf:
                    row[q] = row[p] + 1
                    queue.append(q)
        distances.append(tuple(row))
    embedding = _embedding(internal, adjacency, m)
    if embedding is not None:
        return ({q: nodes[embedding[i]] for i, q in enumerate(logical)},
                [(op[0], *(nodes[embedding[q]] for q in op[1:])) for op in internal])
    best = _fallback(internal, adjacency, distances, m)
    if all(value < inf for row in distances for value in row):
        candidates = [(beam_width, lookahead_weight, False)]
        if portfolio:
            candidates.append((min(beam_width, 500), 0.55, True))
        for width, weight, critical in candidates:
            candidate = _search(internal, adjacency, distances, m, width, weight, critical=critical)
            if (candidate.cost, candidate.swaps) < (best.cost, best.swaps):
                best = candidate
    # A logical state assigned at first use occupied an unused token. Reverse
    # all SWAPs to recover where that token was at circuit input.
    initial = list(best.positions)
    for op in reversed(best.route):
        if op[0] == 'SWAP':
            left, right = op[1:]
            initial = [right if p == left else left if p == right else p for p in initial]
    placement = {q: nodes[initial[i]] for i, q in enumerate(logical)}
    routed = [(op[0], *(nodes[p] for p in op[1:])) for op in best.route]
    return placement, routed
