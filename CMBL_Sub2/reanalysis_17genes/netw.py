import json,sys,networkx as nx
sys.stdout.reconfigure(encoding='utf8')
def metrics(pairs):
    G=nx.Graph(); G.add_edges_from(pairs)
    deg=dict(G.degree()); bt=nx.betweenness_centrality(G,normalized=True)
    cl=nx.closeness_centrality(G,wf_improved=False)
    rows=sorted(G.nodes,key=lambda g:(-deg[g],-bt[g]))
    return G,[(g,deg[g],round(bt[g],4),round(cl[g],4)) for g in rows]
for n in (20,17):
    o=json.load(open('pipeline_%dgenes.json'%n,encoding='utf8'))
    for key in ('consensus','consensus_signcons'):
        if key not in o: continue
        pairs=[tuple(p) for p in o[key]]
        G,rows=metrics(pairs)
        print(f'== {n} genes | {key}: {len(pairs)} pairs, {G.number_of_nodes()} genes, components {nx.number_connected_components(G)}')
        for r in rows: print('   %-8s degree %d  betweenness %.4f  closeness %.4f'%r)
        hubs6=[r[0] for r in rows if r[1]>=6]; hubs5=[r[0] for r in rows if r[1]>=5]
        print('   degree>=6 hubs:',hubs6,'| degree>=5:',hubs5)
