function pat = updatepatchdata(pat, varargin)
% updates patch data if a change is specified
    p = inputParser;
    addParameter(p, 'Areas', [])
    addParameter(p, 'Alphas', [])
    addParameter(p, 'Epsilons', [])
    addParameter(p, 'FaceNormals', [])
    addParameter(p, 'Reference', [])
    addParameter(p, 'Centroids', [])
    parse(p, varargin{:});
    
    areas = p.Results.Areas;
    alphas = p.Results.Alphas;
    epsilons = p.Results.Epsilons;
    norms = p.Results.FaceNormals;
    ref = p.Results.Reference;
    centroids = p.Results.Centroids;

    numfaces = size(pat.Faces, 1);

    % for each face, add associated parameters

    if isempty(pat.UserData) == 1
        pat.UserData = struct('Areas', areas,...
            'Alphas', alphas, 'Epsilons', epsilons, ...
            'FaceNormals', norms, 'Reference',...
            ref, 'Centroids', centroids);
    end

    if not(isempty(areas))
        % sets the areas
        pat.UserData.Areas = areas;
    else
        % sets the areas if they weren't given and not already set
        for i = 1:numfaces
            faces_t = pat.Faces(i, :)';
            pat.UserData.Areas(i) =...
                getarea(pat.Vertices(faces_t(~isnan(faces_t)), :));
        end
    end

    if not(isempty(alphas))
        % sets the alphas
        pat.UserData.Alphas = alphas;
    else
        % sets the alphas if they weren't given and not already set
        warning("Alphas not specified in updatepatchdata after" +...
            " initilization. Default values set to 0.")
        pat.UserData.Alphas = zeros(numfaces, 1);
    end

    if not(isempty(epsilons))
        pat.UserData.Epsilons = epsilons;
    else
    % sets the epsilons if they weren't given and not already set
        warning("Epsilons not specified in updatepatchdata after" +...
            " initilization. Default values set to 0.")
        pat.UserData.Epsilons = zeros(numfaces, 1);
    end
    
    if not(isempty(centroids))
        % sets the centroids
        pat.UserData.Centroids = centroids;
    else
        % sets the centroids if they weren't given and not already set
        pat.UserData.Centroids = zeros(numfaces, 3);
        for i = 1:numfaces
            faces_t = pat.Faces(i, :);
            pat.UserData.Centroids(i,:) =...
                shapecenter(pat.Vertices(faces_t(~isnan(faces_t)), :));
        end
    end

    if not(isempty(ref))
        pat.UserData.Reference = ref;
    else
        % sets the reference if it wasn't given and not already set
        if numfaces == 1
            warning("Patch only has a single face. Reference not" + ...
                " specified in updatepatchdata after initilization." + ...
                " Default reference is a row vector of NaN.")
            pat.UserData.Reference = ref;
        else
            % sets the reference to the centroid of the shape
            pat.UserData.Reference = centroid(pat.Vertices);
        end
    end

    if not(isempty(norms))
        % sets the normal vectors
        pat.UserData.FaceNormals = norms;
    else
        % sets the normal vectors if they weren't given and not already set
        if size(pat.Faces, 1) > 1
            pat.UserData.FaceNormals =...
                calcNormals(pat.Faces, pat.Vertices, pat.UserData.Reference);
        else
            pat.UserData.FaceNormals =...
                calcNormalSimple(pat.Vertices, 2, pat.UserData.Reference);
        end
    end

    % calculates the normal vectors and sets the patch normal vectors
    if isempty(pat.FaceNormals) == 1
        parent_object = pat.Parent;
        parent_type = parent_object.Type;
        iter = 0;
        while strcmp(parent_type, 'axes') == 0  &...
                strcmp(parent_type, 'figure') == 0 & iter < 10
            parent_object = parent_object.Parent;
            parent_type = parent_object.Type;
            iter = iter + 1;
        end

        if and(not(strcmp(parent_type, 'figure') == 1),...
                not(strcmp(parent_type, 'axes') == 1))
            error('Could not reach axes object');
        end
        
        set(pat, 'FaceNormals', pat.UserData.FaceNormals)
    end
    
end