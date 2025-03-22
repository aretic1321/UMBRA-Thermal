function out = valid_d(in, R_s)
    if or(isnumeric(in), isSymType(sym(in), 'number'))
        out = double(in) > R_s;
    elseif isSymType(in, 'variable')
        assume(symvar(in) > R_s)
        out = true;
    else
        out = false;
    end
end