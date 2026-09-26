"""Static placement baseline for comparison with joint placement and routing."""
import random
import networkx as nx
from vendor.starter_kit.scorer import score_summary

def annealed_greedy(program, graph):
    from math import exp
    rng=random.Random(260925)
    logical=sorted({q for op in program for q in op[1:]})
    nodes=sorted(graph.nodes)
    distances=dict(nx.all_pairs_shortest_path_length(graph))
    pairs=[(a,b,.985**i) for i,(kind,*wire) in enumerate(program) if kind=='2Q' for a,b in [wire]]
    def energy(pos):return sum(w*(distances[pos[a]][pos[b]]-1) for a,b,w in pairs)
    best=None
    for trial in range(16):
        current=dict(zip(logical,nodes if trial==0 else rng.sample(nodes,len(logical))))
        value=energy(current)
        lowest,initial=value,dict(current)
        for step in range(1200):
            candidate=dict(current);a=rng.choice(logical);p=rng.choice(nodes)
            for b in logical:
                if candidate[b]==p:
                    candidate[b]=candidate[a];break
            candidate[a]=p
            newvalue=energy(candidate)
            temperature=2*(1-step/1200)+.025
            if newvalue<value or rng.random()<exp((value-newvalue)/temperature):
                current,value=candidate,newvalue
            if value<lowest:lowest,initial=value,dict(current)
        pos=dict(initial);inv={p:q for q,p in pos.items()};route=[]
        for op in program:
            if op[0]=='1Q':route.append(('1Q',pos[op[1]]));continue
            a,b=op[1:];path=nx.shortest_path(graph,pos[a],pos[b])
            for left,right in zip(path[:-2],path[1:-1]):
                la,lb=inv.get(left),inv.get(right);inv[left],inv[right]=lb,la
                if la is not None:pos[la]=right
                if lb is not None:pos[lb]=left
                route.append(('SWAP',left,right))
            route.append(('2Q',pos[a],pos[b]))
        score=score_summary(program,graph,initial,route)['score']
        if best is None or score<best[0]:best=score,initial,route
    return best[1:]

