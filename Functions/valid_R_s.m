function out = valid_R_s(in)
    if or(isnumeric(in), isSymType(sym(in), 'number'))
        out = double(in) > 0;
    elseif isSymType(in, 'variable')
        assume(symvar(in) > 0)
        out = true;
    else
        out = false;
    end
end