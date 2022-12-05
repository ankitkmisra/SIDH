import networkx as nx
import matplotlib.pyplot as plt
from sage.all import *

n = 7
m = 2
p = 2**n * 3**m - 1
F.<i> = GF(p**2, 'i', modulus=x^2+1)

PR.<t> = F[]
def H_p(t):
    m = (p-1)//2
    return sum([binomial(m, i)^2 * t^i for i in range(m+1)])

lambdas = [r[0] for r in H_p(t).roots() if r[0] not in [0,1]]

def subgroup_gen(E, l):
     Ptilde, Qtilde = E.gens()
     P, Q = ((p+1)//l) * Ptilde, ((p+1)//l) * Qtilde
     return P, Q
     

def walk_graph_1(E, S, P_s, Q_s, l, d):
	E_new, P_new, Q_new, S_new = E, P_s, Q_s, S
	
	for i in range(d, 0, -1):
		S_step = l**(i-1) * S_new
		phi = EllipticCurveIsogeny(E_new, [E_new(0), S_step])
		E_new = phi.codomain()
		P_new, Q_new, S_new = phi(P_new), phi(Q_new), phi(S_new)
	
	return E_new, P_new, Q_new


def walk_graph_2(E, S, l, d):
	E_new, S_new = E, S
	
	for i in range(d, 0, -1):
		S_step = l**(i-1) * S_new
		phi = EllipticCurveIsogeny(E_new, [E_new(0), S_step])
		E_new = phi.codomain()
		S_new = phi(S_new)
	
	return E_new

E = EllipticCurve(F, [0, F(-lambdas[0]-1), 0, lambdas[0], 0])
PA, QA = subgroup_gen(E, 2**n)
PB, QB = subgroup_gen(E, 3**m)

print()
print("Predeclared values:")
print("Initial elliptic curve:", E)
print("Alice's basis:", PA, QA)
print("Bob's basis:", PB, QB)
print()

import random
kA = random.randint(0, 2**n - 1)

SA = PA + kA*QA
EA, phiA_PB, phiA_QB = walk_graph_1(E, SA, PB, QB, 2, n)

print()
print("Public key computation:")
print("Alice's secret kA:", kA)
print("Alice's public key:", EA, phiA_PB, phiA_QB)
print()

def oracle(E, R, S, E_j):
	return E_j == walk_graph_2(E, R + kA*S, 2, n).j_invariant()

K = 0
for i in range(n-3):
	kB = random.randint(0, 3**m - 1)
	EB, R, S = walk_graph_1(E, PB + kB*QB, PA, QA, 3, m)
	EAB_j = walk_graph_2(EA, phiA_PB + kB*phiA_QB, 3, m).j_invariant()
	theta = int((1/mod(1 + 2**(n-i-1), 2**n)).sqrt())
	tmp = oracle(EB,  theta*(R - 2**(n-i-1) * K * S), theta*(1 + 2**(n-i-1)) * S, EAB_j)
	K = K + 2**i * (tmp == 0)
	print("computed bits:", bin(K))

for i in range(8):
	K_g = K + i * 2**(n-3)
	if(EA == walk_graph_2(E, PA+K_g*QA, 2, n)):
		print("Alice's secret: ", K_g)

