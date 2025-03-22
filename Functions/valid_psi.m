function out = valid_psi(in)
    if or(isnumeric(in), isSymType(sym(in), 'number'))
        out = and(double(in) > 0, double(in) <= pi);
        if double(in) > pi/2
            % throw warning that function is not implemented yet to work
            % with psi > pi/2
            warning("Function not implemented to work with psi>pi/2"...
                + " at the moement.")
        end
    elseif isSymType(in, 'variable')
        assume([symvar(in) > 0, symvar(in) <= pi])
        out = true;
    else
        out = false;
    end
end