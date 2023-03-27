MAX_DEGREE := 10:
MAX_RECURSION_DEPTH := 10:
VERBOSE := false:

show := proc(msg)
    if VERBOSE then
        print(msg);
    end if;
    end proc:


is_strong_cyclic := proc(A, v, step, allowed_elements)
    local n, B, CV, det, i, j, k, c, eq, sol, solutions, failed, nstep, allowed;
    
    n := RowDimension(A);
    if _npassed = 2 then
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
    
    if _npassed <= 3 then
        allowed := {}
    else
        allowed := allowed_elements
    end if;
    
    solutions := [];
    for i from 1 to n do
        for j from 1 to n do
            if evalb(allowed = {}) or evalb([i, j] in allowed) then
                B := Copy(A);
                prolong_matrix_one(B, c, i, j);
                CV := CV_dv(B, v);
                det := sort(collect(Determinant(CV), x), [x]);
                show(B);
                show(det);
                
                eq := tcoeff(det, x) = 0;
                sol := [solve(eq, c)];
                show(sol);
                #if degree(tcoeff(det, x), c) > 1 then
                    #print(tcoeff(det, x))
                #end if;
                
                for k from 1 to numelems(sol) do
                
                    if sol[k] <> 0 then
                        if eval(det, c=sol[k]) = 0 then
                            # prolongation was found!
                            B := Copy(A);
                            prolong_matrix_one(B, sol[k], i, j);
                            
                            print(B);
                            show("not cyclic!");
                            return false;
                        end if;
                        solutions := [op(solutions), [i, j, sol[k]]];
                    end if;
                end do;
            end if;
        end do;
    end do;
    
    
    failed := false;
    for k from 1 to numelems(solutions) do
        B := Copy(A);
        #print(step, solutions[k]);
        i, j, sol := solutions[k, 1], solutions[k, 2], solutions[k, 3];
        prolong_matrix_one(B, sol, i, j);
        
        if false then
        #if degree(B[i, j], x) > MAX_DEGREE then
            failed := true;
        else
            sol := is_strong_cyclic(B, v, nstep - 1, allowed);
            if sol = false then
                show("not cyclic!");
                return false;
            end if;
            if sol = FAIL then
                failed := true;
            end if;
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