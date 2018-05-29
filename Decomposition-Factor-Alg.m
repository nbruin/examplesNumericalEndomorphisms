/* Decomposition example */
/* Lists all factors in the first example by reconstructing from invariants */

pss := [ ];

ps := [ 43531*x + 4458050224128,
    -43531*x - 139314069504,
    43531*x + 373248000 ];
Append(~pss, ps);

ps := [ -43531*x + 99532800000,
    43531*x - 2280960000,
    -43531*x - 18432000 ];
Append(~pss, ps);

ps := [ -43531*x + 1043677052928,
    43531*x - 39239811072,
    -43531*x + 2217738240 ];
Append(~pss, ps);

Xs := [ ];
for ps in pss do
    g2s := [ Roots(p)[1][1] : p in ps ];
    X := HyperellipticCurveFromG2Invariants(g2s);
    X := ReducedMinimalWeierstrassModel(X);
    f, h := HyperellipticPolynomials(X);
    X := HyperellipticCurve(f / LeadingCoefficient(f));
    Append(~Xs, X);
end for;
print "Factors of the Jacobian:";
print Xs;
