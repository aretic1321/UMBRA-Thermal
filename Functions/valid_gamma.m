function out = valid_gamma(in)
    if isnan(in)
        out = true;
    elseif or(isnumeric(in), isSymType(sym(in), 'number'))
        out = and(double(in) >= -pi, double(in) <= pi);
    elseif isSymType(in, 'variable')
        assume([symvar(in) >= -pi, symvar(in) <= pi])
        out = true;
    else
        out = false;
    end
end