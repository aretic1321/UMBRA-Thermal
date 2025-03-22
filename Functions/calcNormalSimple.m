function [N, order] = calcNormalSimple(vertices, dim, varargin)
% calcNormalSimple finds a normal vector to a set of vertices, where the
% coordinates of each vertex is along dim (1 for row, 2 for column).
% A reference point must be included if you want the normal to face away
% from the point. If there are more than three vertices, it assumes all
% vertices are in the same plane and the shape is convex.
% Also outputs the order of the vertices that would calculate the normals
% anti-clockwise from the perspective of the reference if a reference is
% given. Otherwise, the order is the original order.
    
    p = inputParser;
    addRequired(p, 'vertices')
    addRequired(p, 'dim')
    %{
    if dim == 1
        defref = NaN(3, 1);
    elseif dim == 2
        defref = NaN(1, 3);
    end
    %}
    addOptional(p, 'Reference', [])
    parse(p, vertices, dim, varargin{:});
    
    ref = p.Results.Reference;
    
    mid = shapecenter(vertices, dim);

    % get the vertex coordinates
    if dim == 1
        v1 = vertices(:, 1);
        v2 = vertices(:, 2);
        v3 = vertices(:, 3);
    elseif dim == 2
        v1 = vertices(1, :);
        v2 = vertices(2, :);
        v3 = vertices(3, :);
    end
    % compute two edge vectors
    e1 = v2 - v1;
    e2 = v3 - v2;
    
    % compute the normal vector
    n = cross(e1, e2, dim);
    
    
    if dim == 1
        order = 1:size(vertices, 2);
    else
        order = 1:size(vertices, 1);
    end

    % normalize the normal vector
    if isempty(ref)
        N = n ./ dot(n, n, dim).^(1/2);
    elseif any(isnan(ref))
        N = n ./ dot(n, n, dim) .^(1/2);
    else
        ch = dot(n, mid - ref, dim); % sign gives direction
        if ch < 0
            multip = -1;
            order = fliplr(order);
        else
            multip = 1;
        end
        
        N = n ./ norm(n) .* multip;
    end
end