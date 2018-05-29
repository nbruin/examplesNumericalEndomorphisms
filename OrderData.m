/* Auxiliary function for describing endomorphism rings */

function OrderData(gensHom)
    /* Creation of relevant algebras */
    g := #Rows(gensHom[1]) div 2;
    /* Ambient matrix algebra, plus generators of the endomorphism ring */
    A := Algebra(MatrixRing(Rationals(), 2*g));
    GensA := [ A ! Eltseq(genHom) : genHom in gensHom ];
    /* As a subalgebra */
    B := sub<A | GensA>; GensB := [ B ! gen : gen in GensA ];
    /* As an associative algebra */
    C := AssociativeAlgebra(B); GensC := [ C ! gen : gen in GensB ];

    EndoAlg := [* *]; EndoDesc := [* *];
    EndoAlgQQ, EndoDescQQ := EndomorphismAlgebraQQBase(C);
    Append(~EndoAlg, EndoAlgQQ); Append(~EndoDesc, EndoDescQQ);
    EndoAlgZZ, EndoDescZZ := EndomorphismAlgebraZZBase(C, GensC);
    Append(~EndoAlg, EndoAlgZZ); Append(~EndoDesc, EndoDescZZ);
    EndoAlgRR, EndoDescRR := EndomorphismAlgebraRRBase(C, EndoDescQQ);
    Append(~EndoAlg, EndoAlgRR); Append(~EndoDesc, EndoDescRR);
    return EndoAlg, EndoDesc;
end function;

