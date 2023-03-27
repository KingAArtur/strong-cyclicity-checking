with(LinearAlgebra):
with(SolveTools):

# returns a matrix DV[n x m] = [v | Dv | D2v | ... | D(m-1)v], m = n by default
CV_dv := proc(A, v, m)
    local n, i, j, k, W, AT, width;
    
    n := RowDimension(A);
    if _npassed > 2 then
        width := m
    else
        width := n
    end if;
    
    AT := Transpose(A);
    W := Matrix(n, width);
    for i from 1 to n do
        W[i, 1] := v[i]
    end do;
    
    for j from 2 to width do
        for i from 1 to n do
            W[i, j] := simplify(diff(W[i, j - 1], x) + add(AT[i, k] * W[k, j - 1], k = 1 .. n))
        end do
    end do;
    
    return(W);
    end proc:


# checks if a given vector v is cyclic for system y' = Ay
CV_is_cyclic := proc(A, v)
    local det, answer;
    
    det := Determinant(CV_dv(A, v));
    if det = 0 then
        answer := false
    else
        answer := true
    end if;
    
    return(answer);
    end proc:


# companion matrix for system y' = Ay
CV_companion_matrix := proc(A, v)
    return(MatrixInverse(Transpose(CV_dv(A, v))));
    end proc:


# returns linear scalar differential n-order equation
#   y'(n) = ... for given system z' = Az
# and a matrix B:
#   z = B * [y, y', ..., y'(n-1)]
# using given cyclic vector v
CV_to_scalar_diffur := proc(A, v)
    local n, i, j, DV, C, sys, vars, coefs, eq, y;
    
    n := RowDimension(A);
    DV := CV_dv(A, v, n + 1);
    C := CV_companion_matrix(A, v);
    
    vars := [seq(a[j - 1], j = 1 .. n)];
    sys := [seq(1, i = 1 .. n)];
    for i from 1 to n do
        eq := 0;
        for j from 1 to n do
            eq := eq + DV[i, j] * vars[j]
        end do;
        eq := eq - DV[i, n + 1];
        sys[i] := eq;
    end do;
    
    coefs := Linear(sys, vars);
    
    eq := 0;
    eq := eq + rhs(coefs[1]) * y(x);
    for i from 2 to n do
        eq := eq + rhs(coefs[i]) * diff(y(x), x$(i - 1))
    end do;
    
    return eq = diff(y(x), x$n), C;
    end proc: