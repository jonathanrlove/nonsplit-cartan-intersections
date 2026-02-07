/* Verify claims from accompanying paper */

load "intersection.m";



// Compute intersection numbers of CM points on the provided 
// models of X_{ns}^+(p), and compare each with the prediction of 
// Thereom B when applicable.

// Note that there WILL be disagreement with Theorem B
// for p = 5 and D1 = -3, as there exists a non-Heegner point 
// with CM by this discriminant.

for p in [2,3,5,7,11,13] do
    print "Testing p =",p;
    AllIntersections(p);
end for;



// Check smoothness of models at Fq points

for pair in [<11,2>,<11,3>,<13,2>,<13,3>,<23,2>] do
    p,q := Explode(pair);
    smooth := SmoothAtFqPoints(p,q);
    printf "X_{ns}^+(%o) model smooth at F_{%o}-points: %o\n", p,q,smooth;
end for;
