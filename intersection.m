/* Verify claims from accompanying paper */

load "intersection.m";

print "
       Test 1: Compute intersection numbers of CM points on the provided
       models of X_{ns}^+(p), and compare each with the prediction of
       Thereom B when applicable. 
       
       Note that there WILL be disagreement with Theorem B
       for p = 5 and D1 = -3, as there exists a non-Heegner point
       with CM by this discriminant.
       ";

for p in [2,3,5,7,11,13] do
    print "Testing p =",p;
    AllIntersections(p);
end for;

print "
       Test 2: Find all triples (D1,D2,p) with D1,D2 of class number 1
       and p >= 11 such that S(D1,D2,1) has a multiple of p^2.
       ";

for D1 in ClNum1 do
    for D2 in ClNum1 do
        if D1 gt D2 then
            for M in S(D1*D2, 1) do
                if #[pp[1] : pp in Factorization(M) | pp[1] ge 11 and pp[2] ge 2] ge 1 then
                    printf "%o in S(%o*%o, 1)\n",primefactprint(M), -D1, -D2;
                end if;
            end for;
        end if;
    end for;
end for;



print "
       Test 3: Check that the provided model of X_{ns}^+(p)
       is smooth at F_q-points
       ";

for pair in [<11,2>,<11,3>,<13,2>,<13,3>,<23,2>] do
    p,q := Explode(pair);
    smooth := SmoothAtFqPoints(p,q);
    printf "X_{ns}^+(%o) model smooth at F_{%o}-points: %o\n", p,q,smooth;
end for;
