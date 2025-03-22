function mid = shapecenter(vertices, varargin)
% finds the center of a shape assuming the vertices are coplanar and the
% shape is convex
    p = inputParser;
    addRequired(p, 'vertices')
    addOptional(p, 'dim', 2)
    parse(p, vertices, varargin{:});
    dim = p.Results.dim;
    
    if dim == 1
        oppdim = 2;
    elseif dim == 2
        oppdim = 1;
    end

    % performs triangulation around a mean vertex
    if dim == 1
        vertices = vertices.';
    end
    meanvert = mean(vertices, oppdim);
    C = cat(oppdim, meanvert, vertices);
    X = zeros(size(C, 1)-1, 3);
    for i = 1:(size(C, 1)-2)
        X(i, :) = [1 i+1 i+2];
    end
    X(size(C, 1)-1, :) = [1 size(C, 1) 2];

    % computes the center of the triangles and areas
    cents = zeros(size(X, 1), 3);
    areas = zeros(size(X, 1), 1);
    for i = 1:size(X, 1)
        cents(i, :) = mean(C(X(i, :),:), oppdim);
        areas(i) = getarea(C(X(i, :), :));
    end

    % computes the center point of the shape
    mid = sum(cents.*areas, 1)./sum(areas);
    if dim == 1
        mid = mid.';
    end
end