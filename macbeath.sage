#Example 5.2 (ANTS version of the paper)

#Verification that the curve below is a model for the MacBeath curve.
#Source of model: https://arxiv.org/pdf/1703.01869.pdf 
R.<x,y> = QQ[]
f = 1 + 7*x*y + 21*x^2*y^2 + 35*x^3*y^3 + 28*x^4*y^4 + 2*x^7 + 2*y^7
prec = 100

C = Curve(f)
S = C.riemann_surface(prec = prec)
print "Curve:"
print S
B = S.cohomology_basis(2)
print "Basis of cohomology:"
print B
RG = S.symplectic_isomorphisms()
#we find 1008 automorphisms numerically
assert len(RG) == 1008

#in order to recognize these algebraically we consider their action
#on the tangent space (and hence also on the cotangent space)
Tnum = S.tangent_representation_numerical(RG)

#since our chosen basis of the cotangent space is in terms of differentials
#on the curve that are actually defined over the base field, this
#representation should be acting via transformations with algebraic coefficients.
#since these matrices are a bit large, we use a leap of faith and assume
#they are defined over the cyclotomic field of order 7.
#we can then find probably expressions by finding approximate integer linear dependencies
#between the value and a basis.
k.<zeta> = CyclotomicField(7)
Cr = S._CC
zC = Cr(zeta)
R.<t> = PolynomialRing(QQ)
k.<zeta> = NumberField(t^6 + t^5 + t^4 + t^3 + t^2 + t + 1)
def lincomb(v):
    #w = pari.lindep([v, 1, zC, zC^2, zC^3, zC^4, zC^5])
    w = vector(k, pari.lindep([v, 1, zC, zC^2, zC^3, zC^4, zC^5]).sage())
    if w[0] != 0 :
        return k(list(-w[1:]/w[0]))

#this seems to reconstruct plausible matrices. For instance, one can check they are all
#of finite order
Talg = [ matrix(k, 7,7, [lincomb(a) for a in g.list()]) for g in Tnum]

#in order to get a representation of the automorphisms on the plane 
#model, we do the following. We see that (1:x:y) occur as coordinates 
#in the canonical system specified by our chosen homology basis.

omega = S._differentials
assert omega[4]/omega[0] == x and omega[3]/omega[0] == y

#that means we can compute a rational map P^2->P^2 that should induce the
#appropriate automorphism on C

def mat_to_subs(M):
    w = M*vector(omega)
    return (w[4]/w[0], w[3]/w[0])

#in order to check this we can take the rational map determined by the 
#automorphism and pull back C along it. We can check that C lies in the preimage.

#T = Talg[ ZZ.random_element(len(Talg)) ]
#print "Checking pullback under first element...:"
#w = mat_to_subs(T)
#g = f(*w)
#assert g.numerator() % f == 0
#print "done."

T0 = Talg[0]
print "Checking pullback under first generator..."
w = mat_to_subs(T0)
g = f(*w)
assert g.numerator() % f == 0
print "done."

Tds = [ T for T in Talg if diagonal_matrix(T.diagonal()) == T ]
Tmin = Tds[13]
Td = Tds[0]

print "Checking pullback under second generator..."
w = mat_to_subs(Td)
g = f(*w)
assert g.numerator() % f == 0
print "done."

# We do something stupid because Sage fails to calculate the cardinality
def closure(L):
    print len(L)
    print Tmin in L
    if len(L) > 504:
        return L
    cl = L
    for x1 in L:
        for x2 in L:
            if not x1*x2 in L:
                return closure([ x1*x2 ] + L)
    return L

#print mat_to_subs(Td)
#print mat_to_subs(T0)

# To verify generation claim in paper:
#print Tmin
#cl = closure([ Td, T0 ])
