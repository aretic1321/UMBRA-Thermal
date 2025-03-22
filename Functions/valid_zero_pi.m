function out = valid_zero_pi(in)
    if or(isnumeric(in), isSymType(sym(in), 'number'))
        out = and(double(in) >= 0, double(in) <= pi);
    elseif isSymType(in, 'variable')
        assume([symvar(in) >= 0, symvar(in) <= pi])
        out = true;
    else
        out = false;
    end
end