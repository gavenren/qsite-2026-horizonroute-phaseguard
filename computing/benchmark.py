"""Run public scores, a stronger baseline, ablation, and unseen-case checks."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import platform
import random
import sys
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'vendor'))
import networkx as nx
from starter_kit.hardware import build_hardware_graph, HARDWARE_POSITIONS
from starter_kit.benchmarks import BENCHMARKS
from starter_kit.baseline_routing import solve as official_baseline
from starter_kit.scorer import score_summary
from solver import solve
from baselines import annealed_greedy


def random_greedy(program, graph, trials=100, seed=260925):
    """Best of 100 independent random placements, with shortest-path routing."""
    rng = random.Random(seed)
    logical = sorted({q for op in program for q in op[1:]})
    best = None
    for _ in range(trials):
        initial = dict(zip(logical, rng.sample(list(graph), len(logical))))
        pos = dict(initial)
        inverse = {p: q for q, p in pos.items()}
        routed = []
        for op in program:
            if op[0] == '1Q':
                routed.append(('1Q', pos[op[1]]))
                continue
            a, b = op[1:]
            path = nx.shortest_path(graph, pos[a], pos[b])
            for left, right in zip(path[:-2], path[1:-1]):
                la, lb = inverse.get(left), inverse.get(right)
                inverse[left], inverse[right] = lb, la
                if la is not None:
                    pos[la] = right
                if lb is not None:
                    pos[lb] = left
                routed.append(('SWAP', left, right))
            routed.append(('2Q', pos[a], pos[b]))
        summary = score_summary(program, graph, initial, routed)
        if best is None or summary['score'] < best[0]:
            best = summary['score'], initial, routed
    return best[1:]


def measured(label, function, program, graph):
    start = perf_counter()
    placement, routed = function(program, graph)
    elapsed = perf_counter() - start
    result = score_summary(program, graph, placement, routed)
    if not result['valid']:
        raise AssertionError((label, result['message']))
    result.update(method=label, seconds=round(elapsed, 6), placement=placement, route=routed)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--width', type=int, default=900)
    parser.add_argument('--heldout', type=int, default=4)
    args = parser.parse_args()
    graph = build_hardware_graph()
    rows = []
    cases = {}
    for name, program in BENCHMARKS.items():
        methods = [
            ('Official baseline', official_baseline),
            ('100 random placements', random_greedy),
            ('Static placement search', annealed_greedy),
            ('Small beam (40)', lambda p, g: solve(p, g, beam_width=40, portfolio=False)),
            ('HorizonRoute', lambda p, g: solve(p, g, beam_width=args.width)),
        ]
        cases[name] = {'program': program, 'results': []}
        for method, function in methods:
            result = measured(method, function, program, graph)
            cases[name]['results'].append(result)
            rows.append({'case': name, **{k: result[k] for k in ('method', 'valid', 'swap_count', 'depth', 'score', 'seconds')}})
            print(name, method, result['score'], f"{result['seconds']:.3f}s", flush=True)
    heldout = []
    for seed in range(4101, 4101 + args.heldout):
        rng = random.Random(seed)
        program = [('2Q', *rng.sample(range(8), 2)) for _ in range(20)]
        program.insert(5, ('1Q', 7))
        sample = {'seed': seed, 'program': program, 'results': []}
        for method, function in [('Official baseline', official_baseline), ('100 random placements', random_greedy), ('Static placement search', annealed_greedy),
                                 ('HorizonRoute', lambda p, g: solve(p, g, beam_width=args.width))]:
            result = measured(method, function, program, graph)
            sample['results'].append(result)
            print('heldout', seed, method, result['score'], flush=True)
        heldout.append(sample)
    destination = ROOT / 'results'
    destination.mkdir(exist_ok=True)
    output = {'source_commit': 'a45b04211d4da27dca3c79dbe7009ce58e85dc48',
              'python': platform.python_version(), 'platform': platform.platform(), 'networkx': nx.__version__,
              'beam_width': args.width, 'solver_sha256': hashlib.sha256((ROOT/'solver.py').read_bytes()).hexdigest(),
              'hardware_edges': list(graph.edges), 'hardware_positions': HARDWARE_POSITIONS,
              'public': cases, 'heldout': heldout}
    (destination/'benchmark_results.json').write_text(json.dumps(output, indent=2) + '\n')
    with (destination/'scores.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print('Saved results in', destination, flush=True)

if __name__ == '__main__':
    main()
