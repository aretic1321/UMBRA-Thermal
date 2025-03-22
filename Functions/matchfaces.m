function fir_i_in_sec = matchfaces(fir_faces, fir_vertices, sec_faces, sec_vertices)
% matchfaces matches the first face indicies and puts them in an array
% where each index in the array corrosponds to the index of the second
% face.
    if size(fir_faces,2) > size(sec_faces,2)
        % identifies that a face in the first set has more vertices
        fir_big = true;
        %{
        big_verts = fir_vertices;
        big_faces = fir_faces;
        smal_verts = sec_vertices;
        smal_faces = sec_faces;
        %}
    else
        % identifies that a face in the first set has the same or less
        % amount of vertices
        fir_big = false;
        %{
        big_verts = sec_vertices;
        big_faces = sec_faces;
        smal_verts = fir_vertices;
        smal_faces = fir_faces;
        %}
    end

    if size(fir_faces,1) < size(sec_faces,1)
        less_faces = 1; % less faces in first set
    else
        less_faces = 2; % more or equal num of faces in first set assuming it can't have more
    end
    fir_i_in_sec = zeros(size(sec_faces,1), 0);
    js = 1:size(sec_faces,1);
    for i = 1:size(fir_faces,1) % loops across faces
        for j = js % loops across faces
        %for j = 1:size(sec_faces,1) % loops across faces
            vid1 = fir_faces(i, :)';
            vid2 = sec_faces(j, :)';
            
            if fir_big
                big_verts = fir_vertices(vid1(~isnan(vid1)), :);
                smal_verts = sec_vertices(vid2(~isnan(vid2)), :);
            else
                big_verts = sec_vertices(vid2(~isnan(vid2)), :);
                smal_verts = fir_vertices(vid1(~isnan(vid1)), :);
            end
            % finds if the first and second face have matching vertices or
            % if they exist in the same plane
            if all(ismember(smal_verts, big_verts, "rows") == 1)
                js = js(not(j == js));
                fir_i_in_sec(j) = i;
                break
            elseif less_faces == 1
                % gets vectors between vertices
                vec_chk = big_verts(1:size(smal_verts, 1), :) -...
                    smal_verts;

                % gets the normal vector of the small face
                [N, ~] = calcNormalSimple(smal_verts, 2);
                
                % if the dot product of the normal and vectors between
                % vertices is 0, then the faces are coplanar
                if all(vec_chk*N.' == 0)
                    js = js(not(j == js));
                    fir_i_in_sec(j) = i;
                    break
                end
            end
        end
    end
end