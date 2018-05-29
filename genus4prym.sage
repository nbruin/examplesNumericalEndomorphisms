# Example 5.3 (ANTS version)
#
# Check that a construction by W.P. Milne (1921) that produces a plane
# quartic from a canonical genus 4 curve (with some extra data) indeed
# corresponds to the prym construction.

#We give plane models of the relevant curves
Qyw.<y,w>=QQ[]
C=Curve(y^4*w^2 - y^3*w^3 + y^2*w^4 + 2*y^4*w - y^3*w^2 + 2*y*w^4 + y^4 - 2*y^2*w^2 + y*w^3 + w^4 - y^2*w - y*w^2 + y^2 + 2*y*w + w^2)
Quv.<u,v>=QQ[]
D=Curve(u^4*v^4 - 3*u^4*v^2 + u^4 - u^3*v^3 - 2*u^3*v + u^2*v^2 - u^2 + 3*u*v^3 + 2*u*v + v^4 + v^2 + 1)
Qst.<s,t>=QQ[]
F=Curve(5*s^4 + 28*s^3*t + 28*s^3 + 47*s^2*t^2 + 76*s^2*t + 44*s^2 + 34*s*t^3 + 82*s*t^2 + 66*s*t + 18*s + 16*t^4 + 34*t^3 + 32*t^2 + 18*t + 1)

#Construct the Riemann surfaces
SC=C.riemann_surface(prec=100)
SD=D.riemann_surface(prec=100)
#Singular has several algorithms for computing a canonical basis; they don't
#work well for all input. Here we need "algorithm 2".
SD.cohomology_basis(2)
SF=F.riemann_surface(prec=100)

#some hard work
%time PC=SC.period_matrix()
%time PD=SD.period_matrix()
%time PF=SF.period_matrix()

#Check numerically that the endomorphism rings of Jac(C) and Jac(F)
#are just ZZ (we're finding one generator)
%time EndJC=SC.endomorphism_basis()
print("Rank of End(JC): %o"%len(EndJC))
%time EndJF=SF.endomorphism_basis()
print("Rank of End(JF): %o"%len(EndJF))

#Check numerically that Jac(C) and Jac(F) can be mapped non-constantly
#into Jac(F)
%time HCD=SC.homomorphism_basis(SD)
print("Rank of Hom(JC,JD): %o"%len(HCD))
%time HFD=SF.homomorphism_basis(SD)
print("Rank of Hom(JF,JD): %o"%len(HFD))
