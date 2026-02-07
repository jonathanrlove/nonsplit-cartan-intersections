ZZ := Integers();
QQ := Rationals();
ClNum1 := [-3, -4, -7, -8, -11, -12, -16, -19, -27, -28, -43, -67, -163];

// ==============================================================
//             Models of Non-split Cartan curves
// ==============================================================

// load data for the provided models of X_{ns}^+(p)

load "modular_curve_db.m";
ModCurves := ModularCurveData();



// return the X_{ns}^+(p) model as a scheme over given base ring

function Xnsplus(p : base:=QQ)
	R, Aeqs, jns_0, jns_1 := Explode(ModCurves[p]);
	PP := ProjectiveSpace(ChangeRing(R,base));
    return Scheme(PP, Aeqs);
end function;

// ==============================================================
//               Arithmetic intersection of points
// ==============================================================

// return a list of dehomogenizations with respect to
// every possible component

function EveryCommonAffine(P1,P2)
	d := #P1;
	assert d eq #P2;
	overlap := [i : i in [1..d] | P1[i] ne 0 and P2[i] ne 0];
	affinecoords := [];
	for j in overlap do
		newP1 := [P1[i]/P1[j] : i in [1..d]];
		newP2 := [P2[i]/P2[j] : i in [1..d]];
		Append(~affinecoords, <newP1, newP2>);
	end for;
	return affinecoords;
end function;



// Compute the arithmetic intersection of two points on X_{ns}^+(p).
// Only works if X_{ns}^+(p) is either P^1 or a plane model.

// The result is only correct if for every prime q, the model at q
// is either regular, or becomes regular after blowing up each
// singularity exactly once.
// (this condition holds e.g. for p=11 at q=11.)

function ArithmeticIntersection(P1,P2,p)

	QH := RationalsAsNumberField();
	ZH := Integers(QH);
	d := #P1;

	C := Xnsplus(p);

	affinelist := EveryCommonAffine(P1,P2);
	intprimes := [];
	uncaught := [];
	L := [];

	for ptpair in affinelist do
		// consider each standard affine chart containing the two points
		newP1, newP2 := Explode(ptpair);
		I := &+[(newP1[i]-newP2[i])*ZH : i in [1..d]]; // The congruence ideal! 
		if I eq 0*ZH then
			return [];
		end if;
		newintprimes := [qq : qq in Factorisation(I) | qq[2] gt 0];
		for qq in newintprimes do
			// only record this prime ideal if both points lie in
			// a common affine chart modulo q
			if &and[Valuation(QH!a, qq[1]) ge 0 : a in newP1] and &and[Valuation(QH!a, qq[1]) ge 0 : a in newP2] then
				if not qq[1] in [a[1] : a in intprimes] then
					Append(~intprimes, qq);
					q := qq[1];
					Fq<a>, red := ResidueClassField(q);
					Cq := ChangeRing(C, Fq);
					P1q := Cq ! [red(newP1[i]) : i in [1..d]];
					P2q := Cq ! [red(newP2[i]) : i in [1..d]];
					assert P1q eq P2q;

					// if points intersect at a singularity mod p, subtract one
					// from the intersection multiplicity (to account for blowup)
					if IsSingular(P1q) then 
						if qq[2]-1 ne 0 then Append(~L, <qq[1], qq[2]-1>); end if;
					else
						Append(~L, qq);
					end if;
				else 
					// we've already studied this prime
					assert qq in intprimes;
				end if;
			else
				// we will need to use a different affine chart to study this prime
				Append(~uncaught, qq);
			end if;
		end for;
	end for;

	// check that every "uncaught" prime was eventually found
	assert #[qq : qq in uncaught | not qq[1] in [a[1] : a in intprimes]] eq 0;

	return &*([Norm(qq[1])^qq[2] : qq in L] cat [1]);
end function;


// ==============================================================
//                    Gross-Zagier formula
// ==============================================================

// The epsilon function from Gross-Kohnen-Zagier 
function epsilon(n,D1,D2 : x:=0)

	a,b := SquarefreeFactorisation(n);
	eps := 1;
	assert GCD(D1,D2) eq 1;
	for q in PrimeDivisors(a) do
		if not IsDivisibleBy(D1,q) then eps := eps*KroneckerSymbol(D1,q);
		else eps := eps*KroneckerSymbol(D2,q);
		end if;
	end for;
	return eps;

end function;



// half size of unit group
function w(D)
	if D eq -3 then return 3; end if;
	if D eq -4 then return 2; end if;
	return 1;
end function;



// Compute arithmetic intersections of Heegner divisors. 
function GZFormula(D1,D2 : N := 1)
	prod := 1;
	for t in [0..Ceiling(Sqrt(D1*D2))-1] do
		if IsDivisibleBy(D1*D2-t^2,(4*N)) then
			M  := ZZ!((D1*D2-t^2)/(4*N));
			DM := &*[d^epsilon(ZZ!(M/d),D1,D2) : d in Divisors(M)];
			prod := prod*DM;
		end if;
	end for;
	
	return (ZZ!prod)^(w(D1)*w(D2));
end function;



// ==============================================================
//                  CM points
// ==============================================================

// Reduce flist to a smaller set of polynomials with the 
// same set of roots by eliminating the i-th variable
function EliminateVar(flist, i)
	newlst := [];
	for j in [2..#flist] do
		Append(~newlst, Resultant(flist[1],flist[j],i));
	end for;
	return newlst;
end function;



// plug 1 into the i-th component of the polynomial f, 
// and 0 into each larger component of f
function Dehom(f, i)
	d := Rank(Parent(f));
	newf := Evaluate(f, i, 1);
	for j in [i+1..d] do
		newf := Evaluate(newf, j, 0);
	end for;
	return newf;
end function;



// Compute the rational CM points on X_{ns}^+(p) of discriminant D.
// only works for modular curves that are given as P^1 or a plane curve.

// Note: this can be done much more simply by setting up the j-invariant
// from X_{ns}^+(p) to the j-line as a map of schemes,
// computing the preimage of a desired j-invariant, and searching for 
// rational points in the preimage. However, the method here
// is much faster.
function CMPoints(p, D)
	assert p in [2,3,5,7,11,13];
	assert D in ClNum1;

    R, Aeqs, jns_0, jns_1 := Explode(ModCurves[p]);
	d := Rank(R);
	assert d le 3;
	
	// this function works by searching for rational points in the preimage of 
	// the j-invariant map. fD(x) = x-j_D where j_D is the desired j-invariant.
	fD  := HilbertClassPolynomial(D);
	RQ := PolynomialRing(QQ, d);
	toQ := hom< R -> RQ | [RQ.i : i in [1..d]] >;

	finallist := [];
	for l in [d..1 by -1] do
        // search in affine chart with l-th coordinate equal to 1,
		// only looking for points where all larger coordinates are 0.
		num := Dehom(jns_0,l);
		den := Dehom(jns_1,l);
		dehomAeqs := [Dehom(f,l) : f in Aeqs];
		if den ne 0 then 
			alleqs := dehomAeqs cat [Numerator(Evaluate(fD, num/den))];
			Xiend := [1] cat [0 : i in [l+1..d]];

			// find all possible preimages of j_D under the dehomogenized
			// j-invariant map, working one variable at a time
			Xioptions := [[]];
			for i in [1..l-1] do
				newXioptions := [];
				for Xi in Xioptions do
					eveqs := alleqs;
					for j in [1..i-1] do
						// evaluate polynomials at the components that have already
						// been determined
						eveqs := [Evaluate(toQ(f), j, Xi[j]) : f in eveqs];
					end for;
					reslst := eveqs;
					if i eq 1 and l eq 3 then
						// if polynomials still have two variables, eliminate one
						reslst := EliminateVar(eveqs, 2);
					end if;

					// the result at this point is a univariate polynomial in the 
					// first unknown variable
					roots := Roots(GCD([UnivariatePolynomial(f) : f in reslst]), QQ);
					for root in roots do
						Append(~newXioptions, Append(Xi, root[1]));
					end for;
				end for;
				Xioptions := newXioptions;
			end for;
			if #Xioptions gt 0 then
				for Xi in Xioptions do
					// check that all the candidate points we found do in fact 
					// lie on the curve and have the correct j-invariant
					if &and[Evaluate(f, Xi cat Xiend) eq 0 : f in alleqs] then
						Append(~finallist, Xi cat Xiend);
					end if;
				end for;
			end if;
		end if;
	end for;
	return finallist;
end function;



// ==============================================================
//                     Logistics / Tests
// ==============================================================

// Check that the provided model of X_{ns}^+(p) is smooth at
// all F_q-points.
// (Magma can check for smoothness everywhere, but this
// fails for p=23)
function SmoothAtFqPoints(p,q)
	C := Xnsplus(p : base := GF(q));
	for P in RationalPoints(C) do
		if IsSingular(P) then return false; end if;
	end for;
	return true;
end function;



// Compute all rational CM points on X_{ns}^+(p).
// returns triples <pt, D, i>, where pt is a list of coordinates
// on the model of X_{ns}^+(p) determined in modular_curve_db.m,
// D is the CM discriminant, and i=1,2,... distinguishes between
// multiple points with the same CM order (almost always i=1)
function CMPointList(p)
	PointList := [];
	R<x> := PolynomialRing(QQ);
	for D in ClNum1 do
		foundpts := CMPoints(p, D);

		// swapping the order in one case for purely aesthetic reasons
		if foundpts eq [[-2,1],[-1,1]] then
			foundpts := [[-1,1],[-2,1]];
		end if;

		PointList cat:= [<foundpts[i],D,i> : i in [1..#foundpts]];
	end for;
	return PointList;
end function;



// prints the prime factorization of a rational number.
// use parameter latex to print latex code
function primefactprint(m : latex:=false)
	if m eq 0 then return "0"; end if;
	
	if m lt 0 then sgn := "-";
	else sgn := "";
	end if;

	m := Abs(m);
	if m eq 1 then return sgn cat "1"; end if;
	totstring := "";
	for q in Factorization(Numerator(m)*Denominator(m)) do
		if totstring ne "" then 
			if latex then totstring cat:= "\\,";
			else totstring cat:= "*";
			end if; 
		end if;
		base := IntegerToString(q[1]);
		e := IntegerToString(Valuation(m, q[1]));
		if latex then
			totstring cat:= base cat "^{" cat e cat "}";
		else
			if e eq "1" then totstring cat:= base;
			else totstring cat:= base cat "^" cat e;
			end if;
		end if;
	end for;
	return sgn cat totstring;
end function;



// Tests whether D1, D2 are coprime fundamental discriminants
// and satisfy the (non-split) Heegner hypothesis with respect to p
function GZapplicable(D1,D2,p)
	return GCD(D1,D2) eq 1 
		   and IsFundamentalDiscriminant(D1) and IsFundamentalDiscriminant(D2) 
		   and KroneckerSymbol(D1,p) eq -1 and KroneckerSymbol(D2,p) eq -1;
end function;



// add tag when there are multiple CM points with same discriminant
function tag(t : latex:=false);
	if t eq 1 then 
		return ""; 
	end if;
	if latex then 
		return "_{" cat IntegerToString(t) cat "}";
	end if;
	return "(" cat IntegerToString(t) cat ")";
end function;



// adjust spacing for ease of reading output
function hfill(D)
	if Abs(D) lt 10 then return "  " cat IntegerToString(D); end if;
	if Abs(D) lt 100 then return " " cat IntegerToString(D); end if;
	return IntegerToString(D);
end function;



// Compute all rational CM points on X_{ns}^+(p),
// and all pairwise intersection numbers.
// If D1 and D2 are coprime fundamental Heegner points,
// also compute the intersection number via the Gross-Zagier formula 
// (Theorem B) and check that the numbers agree.

// Use parameter latex to produce a table of intersection numbers.
// The Theorem B check does not run in this case.
procedure AllIntersections(p : latex:=false)

	ptlist := CMPointList(p);

	// set up table headers and first row
	if latex then
		printf "\\begin{array}{|c|";
		for i in [2..#ptlist] do printf "c"; end for;
		printf "|}\n \\hline\n \\Delta_1, \\Delta_2 &  ";
		for j in [2..#ptlist] do
			P2,D2,t2 := Explode(ptlist[j]);
			printf "%o%o ", D2, tag(t2 : latex:=true);
			if j ne #ptlist then printf " & "; end if;
		end for;
		printf " \\\\\n\\hline\n";
	end if;


	for i in [1..#ptlist-1] do
		P1, D1, t1 := Explode(ptlist[i]);
		tag1 := tag(t1 : latex:=latex);

		// first column entry
		if latex then
			printf "%o%o%o", D1, tag1, &cat([" & " : j in [1..i]]);
		end if;

		for j in [i+1..#ptlist] do

			P2, D2, t2 := Explode(ptlist[j]);
			tag2 := tag(t2 : latex:=latex);

			// compute arithmetic intersection
			AIval := ArithmeticIntersection(P1, P2, p);

			if latex then
				// generate table cell entry
				if not GZapplicable(D1,D2,p) then
					printf "\\color{blue} ";
				end if;
				printf " %o ", primefactprint(AIval:latex:=true);
				if j ne #ptlist then printf " & "; end if;
			else
				// print discriminant pair and intersection number
				printf "(%o%o, %o%o): ", hfill(D1),tag1,hfill(D2),tag2;
				printf "%o", primefactprint(AIval);

				// compare with Theorem B prediction
				if GZapplicable(D1,D2,p) then
					GZval := GZFormula(D1,D2 : N:=p^2);
					if GZval eq AIval then
						print " (agrees with Theorem B)";
					else
						print "\n\nDisagrees with Theorem B!";
						print "Theorem B prediction:",primefactprint(GZval);
						print "";
					end if;
				else
					print "";
				end if;
			end if;
		end for;
		// end table row
		if latex then printf " \\\\\n"; end if;
	end for;
	// end table
	if latex then printf "\\hline\n\\end{array}\n"; end if;
end procedure;
