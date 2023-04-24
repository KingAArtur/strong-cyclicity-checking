# Basic functions for checking cyclicity


# returns a matrix DV[n x m] = [v | Dv | D2v | ... | D(m-1)v], m = n by default
CV_dv := proc(A, v, m)
    local n, i, j, k, W, AT, width;
    
    n := LinearAlgebra[RowDimension](A);
    if _npassed > 2 then
        width := m
    else
        width := n
    end if;
    
    AT := LinearAlgebra[Transpose](A);
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
    
    det := LinearAlgebra[Determinant](CV_dv(A, v));
    if det = 0 then
        answer := false
    else
        answer := true
    end if;
    
    return(answer);
    end proc:


# companion matrix for system y' = Ay
CV_companion_matrix := proc(A, v)
    return(LinearAlgebra[MatrixInverse](LinearAlgebra[Transpose](CV_dv(A, v))));
    end proc:


