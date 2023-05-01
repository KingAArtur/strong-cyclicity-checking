# returns a max degree of X in a given polynom (supports O-notation, for example p(x) = x + O(x^5) -> 1)
max_degree := proc(p) :: integer;
    local result :: integer;
    
    if type(p, '`+`') or type(p, polynom) then
        # "polynom" - polynoms with or without "+", such as "x + 3x^3" or "4x^5"
        # "+" - polynoms with O(x^t) term, such as "2x + O(x^3)"
        if has(p, O) then
            # assuming that O(x^t) terms is max degree in the expression
            result := degree(op(select(has, p, O)), x) - 1;
        else
            if p = 0 then
                result := -1
            else
                result := degree(p, x)
            end if;
        end if;
    else
        # only O(x^t) term without any other terms, such as "O(x^4)"
        result := degree(op(p), x) - 1;
    end if;
    
    return result;
end proc:


# removes O(x^t) term from polynom
remove_o_notation := proc(p) :: polynom;
    local q;

    if has(p, O) then
        q := remove(has, p, O);
    else
        q := p;
    end if;

    return q;
end proc:


# returns a prolongation of a given polynom with a given coefficient (symbolic or numeric) and degree
prolong := proc(p :: polynom, coef, degree :: integer) :: polynom;
    local q;
    q := p + coef * x ^ degree;
    
    return q;
end proc:
