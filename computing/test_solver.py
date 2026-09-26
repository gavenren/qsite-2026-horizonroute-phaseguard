"""Correctness tests use the official scorer, not a replacement validator."""
import random
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent/'vendor'))
import networkx as nx
from solver import solve
from starter_kit.scorer import validate_routed_program

class RouterTests(unittest.TestCase):
    def check_route(self, program, graph, **kwargs):
        placement, routed = solve(program, graph, beam_width=24, **kwargs)
        self.assertEqual(validate_routed_program(program, graph, placement, routed), (True, 'ok'))
        self.assertEqual(set(placement), {q for op in program for q in op[1:]})
        return placement, routed

    def test_empty(self):
        self.assertEqual(self.check_route([], nx.Graph()), ({}, []))

    def test_single_qubit_and_labels(self):
        self.check_route([('1Q', 17), ('1Q', 23), ('1Q', 17)], nx.path_graph([101, 303, 909]))

    def test_mixed_order_and_operand_order(self):
        program = [('1Q', 10), ('2Q', 10, 20), ('1Q', 30), ('2Q', 30, 10), ('2Q', 20, 30), ('2Q', 30, 20)]
        self.check_route(program, nx.path_graph(4))

    def test_full_occupancy_star(self):
        self.check_route([('2Q', 0, q) for q in range(1, 9)], nx.path_graph(9))

    def test_disconnected_requires_swaps(self):
        graph = nx.disjoint_union(nx.path_graph(4), nx.path_graph(4))
        program = [('2Q', 0, 1), ('2Q', 4, 5), ('2Q', 0, 2), ('2Q', 4, 6), ('2Q', 0, 3), ('2Q', 4, 7)]
        self.check_route(program, graph)

    def test_no_capacity(self):
        with self.assertRaises(ValueError):
            solve([('2Q', 0, 1)], nx.empty_graph(1))

    def test_impossible_components(self):
        with self.assertRaises(ValueError):
            solve([('2Q', 0, 1), ('2Q', 1, 2)], nx.empty_graph(3))

    def test_bad_operation(self):
        for program in [[('2Q', 1, 1)], [('SWAP', 1, 2)]]:
            with self.assertRaises(ValueError):
                solve(program, nx.path_graph(3))

    def test_determinism(self):
        program = [('2Q', 0, q) for q in range(1, 7)]
        graph = nx.path_graph(8)
        self.assertEqual(self.check_route(program, graph), self.check_route(program, graph))

    def test_random_cases(self):
        for seed in range(64):
            rng = random.Random(9000 + seed)
            n = rng.randint(4, 10)
            graph = nx.path_graph(n) if seed % 2 else nx.cycle_graph(n)
            graph = nx.relabel_nodes(graph, {q: 11 + 7*q for q in graph})
            width = rng.randint(2, n)
            logical = [100 + 3*q for q in range(width)]
            program = [('2Q', *rng.sample(logical, 2)) for _ in range(rng.randint(6, 20))]
            program.insert(3, ('1Q', rng.choice(logical)))
            with self.subTest(seed=seed):
                self.check_route(program, graph)

if __name__ == '__main__':
    unittest.main(verbosity=2)
