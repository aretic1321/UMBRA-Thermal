function q = plotfacenormal(p)
% input: p is a patch object
    
    cent = zeros(size(p.Faces, 1), 3);
    for facen = 1:size(p.Faces, 1)
        % get the indices of the vertices for the face
        v_is = p.Faces(facen, :)';

        % compute the centroid of the shape
        cent(facen, :) = shapecenter(p.Vertices(v_is(~isnan(v_is)), :));
    end
    if not(isempty(p.UserData))
        norms = p.UserData.FaceNormals;
    else
        norms = calcNormals(p.Faces, p.Vertices, p.UserData.Reference);
    end
    hold on
    q = quiver3(cent(:,1), cent(:,2), cent(:,3),...
        norms(:,1), norms(:,2), norms(:,3));
    hold off
end