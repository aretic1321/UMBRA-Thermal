function [faces,swaplist] = CheckOutwardNormals(F, V)
%  Checks normals are outward facing
%  Re-orders face if not
%
%  Input: faces, F (Nx3)
%         vertices, V (Nx3)
%
%  Output: corrected face list so all normals are outwards facing
%
%  Assumptions: Normals are from clockwise orientation of faces. I.e. Right hand  rule
%               Shape is closed and watertight. Run alphaShape if not.
%               No duplicate faces
%
%  Requirements: Function TriangleRayIntersection by Jarek Tuszynski based on "Fast, minimum storage ray-triangle intersection". Tomas Möller and Ben Trumbore. Journal of Graphics Tools,  1997.
%                https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection
%  
%  Modified by Adam Evans March 2025 for improvements

faces=F; % make copy of faces 
swaplist=[]; % list of faces which reordered
% Pick a new origin obviously outside the object
% This will deal with concave/hollow objects where cannot be sure if a point within
% the bounds of the object is actually inside or outside the object
V(:,1) = V(:,1)-1.1*min(V(:,1)); % actually moved object away from the origin
%figure; pch=patch('Faces',F,'Vertices',V); pch.EdgeColor='none'; camlight; pch.FaceColor=[0.8 0.8 0.8]; axis equal vis3d; grid on; hold on; 
for i=1:size(F,1) % for each face
    % find if vector from origin to face is pointing inward or outward at that face
    p = mean(V(F(i,:),:)) ; % point in center of face (assumes face is triangle)
    %hold on; plot3(V(F(i,:),1),V(F(i,:),2),V(F(i,:),3),'r-') % 2 sides of the triangular face
    % plot3([0 p(1)],[0,p(2)],[0,p(3)],'b-') % line segment from origin to centerpoint of face
    % plot3([0 50*p(1)/norm(p)],[0,50*p(2)/norm(p)],[0,50*p(3)/norm(p)],'b-') % ray from origin to centerpoint of face
    % plot3([0 100*p(1)/norm(p)],[0,100*p(2)/norm(p)],[0,100*p(3)/norm(p)],'b-') 
    %[intsect, d, ~, ~, xcoor] = TriangleRayIntersection ([0 0 0], p,V(F(:,1),:), V(F(:,2),:), V(F(:,3),:), 'linetype','segment'); % Not going to use this in case of problems with rounding error
    [intsect, d, ~, ~, xcoor] = TriangleRayIntersection ([0 0 0], p, V(F(:,1),:), V(F(:,2),:), V(F(:,3),:)); % Function by Jarek Tuszynski based on "Fast, minimum storage ray-triangle intersection". Tomas Möller and Ben Trumbore. Journal of Graphics Tools,  1997.
    % intsect is a boolean where 1 if ray passed through face
    % d is the distance in units of |OP| where intersect 
    a=find(intsect) ; % faces which ray OP intersected. One of them should be i. 
    a=setxor(a,i) ; % remove current face i from the set. Ray should penetrate object so a should contain at least the entrance or exit face 
    b=numel(find(d(a)<1)) ; % find how many faces ray intersected before our face of interest 
    % if rem(b,2) then b is even and OP is facing inward when passes through face F(i). Otherwise it is facing out of the object at point p.
    b=rem(b,2); % 0 if OP inwards, 1 if OP outwards
    % find if normal is in same direction as OP 
    n = cross( V(F(i,2),:)-V(F(i,1),:) , V(F(i,3),:)-V(F(i,1),:) );
    n = n/norm(n); % unit normal of face
    % plot3([p(1) p(1)+10*n(1)],[p(2) p(2)+10*n(2)],[p(3) p(3)+10*n(3)],'m-')
    % reorder face if normal is not outward
    c=dot(n,p) ; % dot product of vectors n and p. If >0 n is in same direction as p and if <0 is in opposite direction
    if ( (c>0 && b==0) || (c<0 && b==1) ) % if (same as p and p inward) or (opposite of p and p outward) 
        % reorder face
        swaplist=[swaplist; i];
        faces(i,:)=[faces(i,2),faces(i,1),faces(i,3)];
    end
end
end % function 
