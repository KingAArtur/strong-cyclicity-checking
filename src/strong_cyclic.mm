$include "cyclic.mm"
$include "polynom.mm"


MAX_RECURSION_DEPTH := 10:
VERBOSE := false:
RESULT_MATRIX := Matrix():


show := proc(msg)
    if VERBOSE then
        print(msg);
    end if;
end proc:


find_roots := proc(expr, variable) :: list;
    local roots;
    
    if degree(expr, variable) <= 4 then
        roots := [solve(expr = 0, variable)];
    else
        roots = [];
    end if;
    
    return roots;
end proc:


is_strong_cyclic := proc(A :: Matrix, v :: list, max_recursion_depth :: integer) :: boolean;
    local actual_max_recursion_depth, n, degree_matrix, A_copied, i, j;
    
    if _npassed <= 2 then
        actual_max_recursion_depth := MAX_RECURSION_DEPTH
    else
        actual_max_recursion_depth := max_recursion_depth
    end if;
    
    show("=====");
    show("Starting calculations...");
    show("max_recursion_depth" = actual_max_recursion_depth);
    show("=====");
    
    n := LinearAlgebra[RowDimension](A);
    
    # calculating degrees for each element and removing O-notation
    degree_matrix := Matrix(n);
    A_copied := LinearAlgebra[Copy](A);
    for i from 1 to n do
        for j from 1 to n do
            degree_matrix[i, j] := max_degree(A_copied[i, j]);
            A_copied[i, j] := remove_o_notation(A_copied[i, j]);
        end do;
    end do;
    
    # clearing result matrix
    RESULT_MATRIX := Matrix():
    
    # calling recursive procedure
    return is_strong_cyclic_process(A_copied, v, actual_max_recursion_depth, degree_matrix);
end proc:


is_strong_cyclic_process := proc(A :: Matrix, v :: list, steps_left :: integer, degree_matrix :: Matrix) :: boolean;
    local n, A_copied, degree_matrix_copied, CV, det, i, j, k, c, root, roots, solutions, failed, result;

    show("=====");
    show("depth" = steps_left);
    show(A);
    
    if steps_left = 0 then
        show(FAIL);
        return FAIL;
    end if;
    
    n := LinearAlgebra[RowDimension](A);
    A_copied := LinearAlgebra[Copy](A);
    degree_matrix_copied := LinearAlgebra[Copy](degree_matrix);
    
    # finding all 1-element prolongations that make tr(det) = 0
    solutions := [];
    for i from 1 to n do
        for j from 1 to n do
            A_copied[i, j] := A_copied[i, j] + c * x ^ (degree_matrix_copied[i, j] + 1);
            CV := CV_dv(A_copied, v);
            det := sort(collect(LinearAlgebra[Determinant](CV), x), [x]);
            
            roots := find_roots(tcoeff(det, x), c);
            show(A_copied);
            show(det);
            show(roots);
            
            # checking all found roots and appending them to the list
            for k from 1 to numelems(roots) do
                if roots[k] <> 0 then
                    if eval(det, c=roots[k]) = 0 then
                        # prolongation was found! => not strong cyclic
                        A_copied[i, j] := A[i, j] + roots[k] * x ^ (degree_matrix_copied[i, j] + 1);
                        
                        RESULT_MATRIX := A_copied;
                        show("not cyclic!");
                        return false;
                    end if;
                    solutions := [op(solutions), [i, j, roots[k]]];
                end if;
            end do;
            
            # discarding changes for this element
            A_copied[i, j] := A[i, j];
        end do;
    end do;
    
    # recursively applying to all collected prolongations
    failed := false;
    for k from 1 to numelems(solutions) do
        i, j, root := solutions[k, 1], solutions[k, 2], solutions[k, 3];
        
        A_copied[i, j] := A_copied[i, j] + root * x ^ (degree_matrix_copied[i, j] + 1);
        degree_matrix_copied[i, j] := degree_matrix_copied[i, j] + 1;
        
        result := is_strong_cyclic_process(A_copied, v, steps_left - 1, degree_matrix_copied);
        degree_matrix_copied[i, j] := degree_matrix[i, j];
        A_copied[i, j] := A[i, j];
        
        if result = false then
            show("not cyclic!");
            return false;
        end if;
        if result = FAIL then
            failed := true;
        end if;
    end do;
    
    if failed then
        show("FAIL");
        return FAIL;
    else
        show("yes, strong cyclic!");
        return true;
    end if;
end proc: