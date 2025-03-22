% CENTROID - computes the center of gravity of a convex polyhedron in any
%            number of dimensions
%
% USAGE: C = centroid(P)
% 
% P = Matrix of convex polyhedron vertices: each row is a vertex, and each
%     column is a dimension.
% C = Row vector of centroid coordinates. Each column is a dimension.
%
% Notes: (1) This function computes the centroid by partitioning into
%            simplices and determining the weighted sum of their centroids.
%        (2) Written in response to a posting on comp.soft-sys.matlab
%
% Michael Kleder, Sep 2005
%
% Modified by Adam Evans, Dec 2024 to work with current version of MATLAB
% 

function C = centroid(P)
    k=convhulln(P);
    if length(unique(k(:)))<size(P,1)
        error('Polyhedron is not convex.');
    end
    if size(P, 2) > 3
        T = delaunayn(P);
    else
        T = delaunay(P);
    end
    n = size(T,1);
    W = zeros(n,1);
    C=0;
    for m = 1:n
        sp = P(T(m,:),:);
        [~,W(m)]=convhulln(sp);
        C = C + W(m) * mean(sp);
    end
    C=C./sum(W);
end