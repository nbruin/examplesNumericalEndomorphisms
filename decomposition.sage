#Example 5.1 (ANTS version of the paper)

#An example of how to use numerical computation of idempotents in endomorphism rings to
#discover non-galois subcovers.

#We consider the following curve:

A.<x,y>=AffineSpace(QQ,2)
C=Curve(4*x^6-54*x^5*y-729*x^4+108*x^3*y^3+39366*x^2-54*x*y^5-531441)

#It has an endomorphism ring (at least numerically) that looks like Z x Z x Z

S=C.riemann_surface(prec=100)
print "Period matrix:"
print S.period_matrix()
E=S.endomorphism_basis()
assert len(E) == 3

#The automorphism group is of order 2 (this is a non-hyperelliptic curve so the automorphism group is
#the quotient by <-1> (and -1 is in the list of length 4 found below)

assert len(S.symplectic_isomorphisms(hom_basis=E)) == 4

#The non-trivial automorphism is easily found by inspection or by looking at the action on the space of
#differentials: it's (x,y) :-> (-x,-y).

#There are only finitely many idempotents, so we determine them by taking a general member of the
#endomorphism ring M and insist that M^2-M=0, and take the coefficient-wise equations. This forms a 0-dimensional
#ideal and we determine the elements.

P.<a,b,c>=QQ[]
M=sum(a*m for a,m in zip(P.gens(),E))
idempotent_combinations=[tuple(v[a] for a in P.gens()) for v in P.ideal((M^2-M).list()).variety()]
idempotents=[sum(a*m for a,m in zip(v,E)) for v in idempotent_combinations]
idems_rank8=[M for M in idempotents if M.rank() == 8]

#We find that there are 3 idempotents with a 2-dimensional kernel in their representation on the cotangent space
#We determine the bases of these kernels.

#Note that we have a representation of a canonical system of this curve 
#in the form of our basis of regular differentials. We can determine the 
#subspaces that correspond to them. If C admits a cover phi: C->D, then 
#phi^*H^0(D,Omega^1_D) will be a kernel of some idempotent in the 
#endomorphism algebra.

pis=S.tangent_representation_algebraic(idempotents)
W=[p.kernel().basis_matrix()*vector(S._differentials) for p in pis if p.rank() == 4]

#We find that these three 2-dimensional spaces span the whole space, so it 
#looks like Jac(C) decomposes into three abelian surfaces.

#Given W subset H^0(C,Omega^1), we can consider the projection 
#C->Proj(W). If W corresponds to phi^*H^0(D,Omega^1) for some subcover, 
#then we will find that this map to Proj(W) factors through phi. We determine the covers.

#Note that in these cases, Proj(W)=P^1, so we just get D expressed as a cover of the projective line
#in various ways. We get the following functions:

print "Functions to PP^1:"
print [w[0]/w[1] for w in W]

#from which we find: [y/x, (y^2 - 81)/(x^2 - 81), y/x]. We determine plane models of D, where one of the coordinates
#is the function t suggested here. From

Qtxy=QQ["t,x,y"]; t = Qtxy('t'); y = Qtxy("y")
print "Resultant:"
print (x*t-y).resultant(Qtxy(S.f),Qtxy('y'))

#we find: 54*t^5*x^6 - 108*t^3*x^6 + 54*t*x^6 - 4*x^6 + 729*x^4 - 39366*x^2 + 531441=0, which indeed describes a genus 6 curve again.

#We can get our hands on the Galois theory of this extension of the t-line using monodromy computations:

G=SymmetricGroup([0..5]).subgroup(Curve(QQ['t,x']((x*t-y).resultant(Qtxy(S.f),y))).riemann_surface().monodromy_group())
assert G.is_transitive()
H=G.stabilizer(0)
assert len(H) == 2 and len(G) == 12
print "Cardinalities of subgroups:"
print [len(g) for g in G.subgroups() if H.is_subgroup(g)]

#gives [2,4,6,12]. So Galois theory tells us that this degree 6 cover of 
#the t-line has a degree 2 and a degree 3 subcover. Indeed, the element

#(5832*t^5 - 11664*t^3 + 5832*t - 432)*x^5 + (157464*t^5 - 314928*t^3 + 
#    157464*t + 67068)*x^3 - 2125764*x

#has minimal polynomial:

#x^2 + 244019119519584*t^5 - 488038239039168*t^3 + 244019119519584*t - 18075490334784

#Describing a genus 2 hyperelliptic curve, admitting the model:

#y^2 = -16*x^5 - 40*x^4 + 32*x^3 + 88*x^2 - 32*x - 23

#The computation above shows that D is a cubic cover of this genus 2 curve. Note that this cover is NOT galois,
#so this curve is NOT obtained by taking a group quotient of D.

#Similarly, the element

#(8503056*t^5 - 17006112*t^3 + 8503056*t - 629856)/(t^2 - 1)*x^4 + 
#    (688747536*t^5 - 1377495072*t^3 + 688747536*t + 63772920)/(t^2 - 1)*x^2 - 
#    1033121304/(t^2 - 1)

#has minimal polynomial

#x^3 - 172909019862142987392*t*x + 3215447857320103140639814785024*t^4 - 3215447857320103140639814785024*t^2

#which is a model for the genus 2 curve:

#y^2 = x^6 - x^5 + 1

#The curve D is a double cover of this curve, so this is a quotient of D.

#The final factors is found by reconstructing from Igusa invariants; see the Magma files in this repository.

#We can now use sage to determine the homomorphisms between the respective complex tori.

QX.<X>=QQ[]
C1=HyperellipticCurve(X^6-X^5+1)
C2=HyperellipticCurve(-16*X^5 - 40*X^4 + 32*X^3 + 88*X^2 - 32*X - 23);
C3=HyperellipticCurve(X^6+3*X^4+3*X^2-X+1);
SC1=C1.riemann_surface(prec=100)
SC2=C2.riemann_surface(prec=100)
SC3=C3.riemann_surface(prec=100)

A1=SC1.homomorphism_basis(S)
Ar1=S.homomorphism_basis(SC1)
A2=SC2.homomorphism_basis(S)
Ar2=S.homomorphism_basis(SC2)
A3=SC3.homomorphism_basis(S)
Ar3=S.homomorphism_basis(SC3)

#A side-effect of LLL is that we should be finding minimal degree maps.
#It's not guaranteed that by asking for isomorphisms we actually find
#duals, but in this case the composition indeed give scalar multiplication
#by a positive integer.

assert Ar1[0]*A1[0] == 2
assert Ar2[0]*A2[0] == 3
assert Ar3[0]*A3[0] == 6

assert Ar1[0].row_space().intersection(Ar2[0].row_space()).dimension() == 0
assert Ar1[0].row_space().intersection(Ar3[0].row_space()).dimension() == 0
assert Ar2[0].row_space().intersection(Ar3[0].row_space()).dimension() == 0

#Together this confirms S is isogenous to the sum of SC1,SC2,SC3
#We can verify that independently:

SUM=SC1+SC2+SC3
ES=SUM.homomorphism_basis(S)
assert sum(ES).det() != 0
