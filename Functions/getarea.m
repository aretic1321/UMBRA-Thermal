function A = getarea(vertices, varargin)
% getarea computes the area of a polygon assuming the vertices are listed
% along the rows unless otherwise specific. dim specifies the dimension
% that the coordinates of a vertex are listed across

    p = inputParser;
    addRequired(p, 'vertices')
    addOptional(p, 'dim', 2)
    parse(p, vertices, varargin{:});
    
    dim = p.Results.dim;

    if dim == 2
        v = vertices - mean(vertices, 1);
        A = area3D(v(:, 1), v(:, 2), v(:, 3));
    elseif dim == 1
        v = vertices - mean(vertices, 1);
        A = area3D(v(1, :), v(2, :), v(3, :));
    end
end
