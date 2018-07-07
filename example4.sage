A.<x,y>=AffineSpace(QQ,2)
C=Curve(3*x^5 + 2*x^3*y^2 - x*y^4 - 9*x^4 + 6*x^2*y^2 - y^4 + 4)

%time S=C.riemann_surface(prec=100)
%time M=S.period_matrix()

print "Period matrix:"
print S.period_matrix()

#we compute the symplectic isomorphisms of the Jacobian
#with its action on the homology:
auts_homology=S.symplectic_isomorphisms()

#and its action on the tangest space
auts_tangent=S.tangent_representation_algebraic(auts_homology)

#Note that (1:x:y) can be consider a sublinear system of the canonical
#system on this curve, and that we can find these plane coordinates in
#our basis for the differentials:

assert [S._differentials[i] for i in [0,1,3]] ==[1,x,y]

#this means we can compute the action of the automorphisms on the planar model
#quite easily:
v=vector(S._differentials)
w=[v*a.transpose() for a in auts_tangent]
auts_planar=[[c[1]/c[0],c[3]/c[0]] for c in w]

print "automorphisms acting on the plane model of the curve:"
print auts_planar

#note that they all act linearly on the affine model!
#we check the automorphisms really do leave the equation invariant
assert all(C.defining_polynomial()(a) == C.defining_polynomial() for a in auts_planar)       

#we select the involutions from the automorphisms. Note that if a is one, then
# so is -a. We only select the non-trivial ones, and we pick the one of larger rank
# (which is 8 in this case)

involutions=[a for a in auts_homology if a^2==1 and (1-a).rank() == 8]

#note that the automorphism group of the curve is the quotient by -1
#so the following pins the group as Sym(3).
assert len(auts_homology)==12 and len(involutions)==3

#each involution leaves a 4-dimensional subspace of the homology fixed.
#the sum of these spaces is 8-dimensional:
V=sum([(1-i).right_kernel() for i in involutions])

#we consider a hyperelliptic curve:
D1=Curve(y^2-(x^6-x^5+1))
SD1=D1.riemann_surface(prec=100)
%time _=SD1.period_matrix()

#and compute the homomorphisms between Jac(D1) and Jac(C)
H1basis=SD1.homomorphism_basis(S)

#the subspace spanned by the images of these homomorphisms agrees with V
#this confirms (numerically) that Jac(D1)^2 is an isogeny factor of Jac(C).
#of course, a groebner-basis type computation can compute quotients of 
#C by these involutions and find that these are isomorphic to D1.

W=sum([h.column_space() for h in H1basis])
assert V==W

#we determine the order 3 elements in the group
order3=[a for a in auts_homology if a^3==1 and a != 1]

assert len(order3)==2 and all( (1-a).rank() == 8 for a in order3)

#We consider another genus 2 curve:

D2=Curve(y^2-(2*x^5+27*x^4-54*x^2+27))
SD2=D2.riemann_surface(prec=100)
%time _=SD2.period_matrix()

#and compute homomorphisms from Jac(D2) to Jac(C) (only 1-dimensional):

H2basis=SD2.homomorphism_basis(S)
assert len(H2basis)==1

#and its image agrees with the fixed space of the order 3 subgroup.
assert H2basis[0].column_space() == (1-order3[0]).right_kernel()

#indeed, with other means one can check that the quotient of C by the
#order 3 automorphism group gives a curve isomorphic to D2.
