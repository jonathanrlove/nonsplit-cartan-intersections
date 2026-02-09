# Non-split Cartan Intersections
Accompanies the paper "Arithmetic intersections on non-split Cartan modular curves" by Jonathan Love, Elie Studnia, and Jan Vonk (arxiv identifier to be added once uploaded).

## Quick start

To verify the claims in section 4 of the accompanying paper, simply run main.m. This runs the following tests.

For each $p \in \\{2, 3, 5, 7, 11, 13\\}$: 
- Take the provided integral model of $X_{\mathrm{ns}}^+(p)$.
- Find all rational CM points on this model.
- Compute all pairwise intersections of these points. The code will print `(D1,D2):m` if points with CM discriminant `D1` and `D2` have intersection number `m`. If there are multiple CM points with the same discriminant, tags `(2)`, `(3)`, etc will be appended to the discriminant to distinguish each point.
- When the CM points have coprime fundamental discriminants that satisfy the Heegner hypothesis ($p$ inert in the CM order), the computed intersection is compared against the value predicted by Theorem B.
- **IMPORTANT NOTE**: there exist two rational points on $X_{\mathrm{ns}}^+(5)$ with CM discriminant $-3$. One of these is not a Heegner point, and the Theorem B checks WILL fail for this point. This is to be expected. The tests should succeed for the Heegner point of the same discriminant.

For each $p \in \\{11, 13, 23\\}$:
- Take an integral model of $X_{\mathrm{ns}}^+(p)$.
- Check that all $\mathbb{F}_q$-points of this model are smooth, for $q=2,3$ (for $p=23$ we only check $q=2$).

## Contents of repository

There are three Magma files. 

### `main.m` 
This file runs the checks described in the previous section.

### `modular_curve_db.m` 
This file contains models for non-split Cartan modular curves $X_{\mathrm{ns}}^+(p)$ for $p \in \\{2, 3, 5, 7, 11, 13, 23\\}$. These models are cut out by a set `model_0` of $n-1$ polynomial equations in $\mathbb{P}^n$, and come equipped with a $j$-invariant map `[map_0_coord_0, map_1_coord_1]` to $\mathbb{P}^1$. By imposing the condition that the $j$-invariant map must be flat over $\mathbb{P}^1_{\mathbb{Z}}$, we obtain a uniquely determined integral model of $X_{\mathrm{ns}}^+(p)$, even in the case that the generic fibre has genus $0$.

The data for these modular curves was obtained from the LMFDB, namely the following pages (last accessed February 7, 2026):

- https://beta.lmfdb.org/ModularCurve/Q/1.1.0.a.1/
- https://beta.lmfdb.org/ModularCurve/Q/3.3.0.a.1/
- https://beta.lmfdb.org/ModularCurve/Q/5.10.0.a.1/
- https://beta.lmfdb.org/ModularCurve/Q/7.21.0.a.1/
- https://beta.lmfdb.org/ModularCurve/Q/11.55.1.b.1/
- https://beta.lmfdb.org/ModularCurve/Q/13.78.3.a.1/
- https://beta.lmfdb.org/ModularCurve/Q/23.253.13.a.1/

For $p\in\\{2,3,5,7\\}$, the values of $p$ such that $X_{\mathrm{ns}}^+(p)\simeq \mathbb{P}^1_{\mathbb{Q}}$, the $j$-map given by the LMFDB was precomposed with a rational automorphism of $\mathbb{P}^1$ to determine a different integral structure.

### `intersection.m`
This file contains methods for computing rational CM points on modular curves and computing arithmetic intersections of these points. See the file for documentation on how to use method.

The main methods are the following:
- `CMPointList(p)` returns a list of all rational CM points on the provided model of $X_{\mathrm{ns}}^+(p)$ (i.e. all rational points for which the $j$-invariant map sends the point to the $j$-invariant of a CM discriminant).
- `ArithmeticIntersection(P1,P2,p)` takes points $P_1,P_2\in X_{\mathrm{ns}}^+(p)(\mathbb{Q})$ as input (provided as a list of coordinates), and returns the product over all primes of $q^{m_q}$, where $m_q$ is the intersection multiplicity of $P_1$ and $P_2$ at $q$.
- `GZFormula(D1, D2 : N:=p^2)` returns the intersection number of Heegner points on $X_{\mathrm{ns}}^+(p)$ with discriminants $D_1$ and $D_2$, computed using Theorem B.
- `AllIntersections(p)` computes all rational CM points on $X_{\mathrm{ns}}^+(p)$ using `CMPointList(p)`, computes all pairwise intersections using `ArithmeticIntersection(P1,P2,p)`, and compares the result to the Theorem B prediction using `GZFormula(D1, D2 : N:=p^2)` whenever the CM discriminants are coprime, fundamental, and satisfy the Heegner hypothesis ($p$ inert). Using the `latex` parameter, this code can also produce tables of intersection numbers.
- `SmoothAtFqPoints(p,q)` checks that the provided integral model of $X_{\mathrm{ns}}^+(p)$ is smooth at all $\mathbb{F}_q$-points.

