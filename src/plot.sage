import networkx as nx
import matplotlib.pyplot as plt
from sage.all import *

p = 2**4 * 3**3 - 1
F.<i> = GF(p**2, 'i', modulus=x^2+1)

def plot(G, clist, label = True):
	nx.draw_spring(G, with_labels=label, node_color=clist)
	plt.show()

def isogeny_graph(N):
	G = nx.Graph()
	E_0 = EllipticCurve(F, [0, 208*i+161, 0, 1, 0])

	visited = {}
	to_visit = [E_0]
	visited[E_0.j_invariant()] = len(visited)

	while(len(to_visit) != 0):
		e = to_visit[0]
		assert(e.is_supersingular())

		for ker_phi_point in e(0).division_points(N):
			if(ker_phi_point.is_zero()):
				continue
			
			phi = EllipticCurveIsogeny(e, [e(0), ker_phi_point])

			G.add_edge(e.j_invariant(), phi.codomain().j_invariant())
			
			if(phi.codomain().j_invariant() not in visited.keys()):
				to_visit = to_visit + [phi.codomain()]
				visited[phi.codomain().j_invariant()] = len(visited)
		
		to_visit = to_visit[1:]

	return G

def walk_graph(E, S, P_s, Q_s, l, d):
	colorset = set()
	E_new, P_new, Q_new, S_new = E, P_s, Q_s, S
	
	for i in range(d, 0, -1):
		S_step = l**(i-1) * S_new
		phi = EllipticCurveIsogeny(E_new, [E_new(0), S_step])
		E_new = phi.codomain()
		P_new, Q_new, S_new = phi(P_new), phi(Q_new), phi(S_new)
		colorset.add(E_new.j_invariant())
	
	return E_new, P_new, Q_new, colorset


G_2 = isogeny_graph(2)
G_3 = isogeny_graph(3)

E = EllipticCurve(F, [0, 329*i+423, 0, 1, 0])
PA, QA = E(100*i + 248, 304*i + 199), E(426*i + 394, 51*i + 79)
PB, QB = E(358*i + 275, 410*i + 104), E(20*i + 185, 281*i + 239)

print()
print("Predeclared values:")
print("Initial elliptic curve:", E)
print("Alice's basis:", PA, QA)
print("Bob's basis:", PB, QB)
print()

clist_2 = ['grey' if node == E.j_invariant() else 'cyan' for node in G_2.nodes()]
clist_3 = ['grey' if node == E.j_invariant() else 'cyan' for node in G_3.nodes()]
plot(G_2, clist_2)
plot(G_3, clist_3)

import random
kA = 12 #random.randint(0, 2**4 - 1)
kB = 3 #random.randint(0, 3**3 - 1)

SA = PA + kA*QA
SB = PB + kB*QB
EA, phiA_PB, phiA_QB, csetA = walk_graph(E, SA, PB, QB, 2, 4)
EB, phiB_PA, phiB_QA, csetB = walk_graph(E, SB, PA, QA, 3, 3)

print()
print("Public key computation:")
print("Alice's secret kA:", kA)
print("Bob's secret kB:", kB)
print("Alice's public key:", EA, phiA_PB, phiA_QB)
print("Bob's public key:", EB, phiB_PA, phiB_QA)
print()

clist_2 = ['orange' if node in csetA else 'lightgreen' if node in csetB else clist_2[i] for i,node in enumerate(G_2.nodes())]
clist_3 = ['orange' if node in csetA else 'lightgreen' if node in csetB else clist_3[i] for i,node in enumerate(G_3.nodes())]
plot(G_2, clist_2)
plot(G_3, clist_3)

SA = phiB_PA + kA*phiB_QA
SB = phiA_PB + kB*phiA_QB
EA_shared, _, _, csetA = walk_graph(EB, SA, phiB_PA, phiB_QA, 2, 4)
EB_shared, _, _, csetB = walk_graph(EA, SB, phiA_PB, phiA_QB, 3, 3)

print()
print("Shared secret key computation:")
print("Alice's shared secret key:", EA_shared.j_invariant())
print("Bob's shared secret key:", EB_shared.j_invariant())
print()

clist_2 = ['mediumpurple' if node in csetA else clist_2[i] for i,node in enumerate(G_2.nodes())]
clist_3 = ['mediumpurple' if node in csetB else clist_3[i] for i,node in enumerate(G_3.nodes())]
plot(G_2, clist_2)
plot(G_3, clist_3)

