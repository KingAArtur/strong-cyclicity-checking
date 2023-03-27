# returns a prolongation of a given polynom (supports O-notation, for example q(x) = x + O(x^5))
prolong := proc(q, c, need_indexes, i, j)
    local t, p, cc;
    p := q;
    
    if type(p, '`+`') or type(p, polynom) then
        if has(p, O) then
            t := degree(op(select(has, p, O)), x);
            p := remove(has, p, O);
        else
            if p = 0 then
                t := 0
            else
                t := degree(p, x) + 1
            end if;
        end if;
    else
        t := degree(op(p), x);
        p := 0;
    end if;
    
    if _npassed <= 2 then
        cc := c
    else
        if _npassed >= 5 then
            cc := c[i, j, t]
        else
            cc := c[t]
        end if;
    end if;
    
    
    p := p + cc * x^t;
    return p;
    end proc:


# prolong one element of a given matrix
prolong_matrix_one := proc(A, c, i, j)
    A[i, j] := prolong(A[i, j], c);
    end proc: