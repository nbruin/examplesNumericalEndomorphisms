# Numerical endomorphism ring calculation for the curve

A.<x,y> = AffineSpace(QQ, 2)
C = Curve(4*x^6 - 54*x^5*y - 729*x^4 + 108*x^3*y^3 + 39366*x^2 - 54*x*y^5 - 531441)
S = C.riemann_surface(prec = 100)

Rs = S.endomorphism_basis()
print Rs
magma.load("OrderData.m")
print magma.OrderData(Rs, nvals = 2)
print "The entry [RR, RR, RR] indicates that we have an order in RR x RR x RR and hence some subring of the maximal order ZZ x ZZ x ZZ. The entry [6, -1] indicates that the index in this maximal order is 6."
