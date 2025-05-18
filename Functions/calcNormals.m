function [normals, order] = calcNormals(faces, vertices, varargin)
% calcNormals finds normal vectors of a face or faces. To find the normal
% vector of a single face, just make sure faces has some array where the
% size of the rows is 1.
    
    p = inputParser;
    addRequired(p, 'Faces')
    addRequired(p, 'Vertices')
    addOptional(p, 'Reference', [])
    
    parse(p, faces, vertices, varargin{:})
    
    ref = p.Results.Reference;
    if size(faces, 1) == 1
        [normals, order] = calcNormalSimple(vertices, 2, ref);
    else
        DT = delaunayTriangulation(vertices);
        [T, X] = freeBoundary(DT);
        TR = triangulation(T, X);
        Ns = faceNormal(TR);
        
        % uses matchFaces to determine which normal vecotrs match to
        % which face, then applies them to the associated original face
        normals = zeros(size(faces,1), 3);
        oldtonew = matchfaces(faces,vertices, T, X);
        order = oldtonew;
        for i = 1:size(faces, 1)
            n_i = find(oldtonew == i, 1, 'first');
            normals(i, :) = Ns(n_i, :);
        end
    end
    % enforces the normals to have a length of 1
    normals = normals./dot(normals, normals, 1).^(1/2);
end