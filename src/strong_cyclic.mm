$include "cyclic.mm"
$include "prolong.mm"

MAX_RECURSION_DEPTH := 10:
VERBOSE := false:

show := proc(msg)
    if VERBOSE then
        print(msg);
    end if;
end proc:


is_strong_cyclic := proc(A :: Matrix, v :: list, degree :: integer, step :: integer) :: boolean;
    local n, B, CV, det, i, j, k, c, eq, sol, solutions, failed, nstep;
    
    if _npassed = 3 then
        nstep := MAX_RECURSION_DEPTH
    else
        nstep := step
    end if;
    
    show("=====");
    show("depth" = nstep);
    show(A);
    
    if nstep = 0 then
        show(FAIL);
        return FAIL;
    end if;
    
    n := LinearAlgebra[RowDimension](A);
    
    # finding all 1-element prolongations that make tr(det) = 0
    solutions := [];
    for i from 1 to n do
        for j from 1 to n do
            B := LinearAlgebra[Copy](A);
            prolong_matrix_one(B, c, i, j, degree + 1);
            CV := CV_dv(B, v);
            det := sort(collect(LinearAlgebra[Determinant](CV), x), [x]);
            show(B);
            show(det);
            
            eq := tcoeff(det, x) = 0;
            sol := [solve(eq, c)];
            show(sol);
            
            for k from 1 to numelems(sol) do
                if sol[k] <> 0 then
                    if eval(det, c=sol[k]) = 0 then
                        # prolongation was found!
                        B := LinearAlgebra[Copy](A);
                        prolong_matrix_one(B, sol[k], i, j, degree + 1);
                        
                        print(B);
                        show("not cyclic!");
                        return false;
                    end if;
                    solutions := [op(solutions), [i, j, sol[k]]];
                end if;
            end do;
        end do;
    end do;
    
    # recursively applying to all collected prolongations
    failed := false;
    for k from 1 to numelems(solutions) do
        B := LinearAlgebra[Copy](A);
        i, j, sol := solutions[k, 1], solutions[k, 2], solutions[k, 3];
        prolong_matrix_one(B, sol, i, j, degree + 1);
        
        sol := is_strong_cyclic(B, v, degree + 1, nstep - 1);
        if sol = false then
            show("not cyclic!");
            return false;
        end if;
        if sol = FAIL then
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