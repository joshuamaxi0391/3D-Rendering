% HW1_Joshua_Maxi.m

%  INITIALIZED VARIABLES TO CHANGE HOW THE 3D OBJECT IS VIEWED AND DISPLAYED IN 2D WORLD SPACE
%─────────────────────────────────────────────────────────────────────────────────────────────

% Object translation, scaling, and rotation (in degrees)
objectRotationX_degree = 45;    
objectRotationY_degree = 0;    
objectRotationZ_degree = 110;     
Sx = 2; Sy = 2; Sz = 2;    
Tx = 10; Ty = -5; Tz = 7; 

% Camera x, y, z positions
cameraPositionX = -20;       
cameraPositionY = -20;      
cameraPositionZ = -10;      

% Have camera look at the origin
camLookAt = [0, 0, 0]; 

% Light position and material diffuse color
lightPosition = [-8, -8, 2]; 
M_diff = [1, 0, 0];

% Perspective variables
fov = 60;    
aspect = 1;  
nearZ = 0.1;   
farZ = 100.0; 

% Filename to determine the 3D object to be rendered
filename = 'mew_lp.raw';

% WORLD MATRIX INITIALIZATION AND CALCULATION
%─────────────────────────────────────────────────────────────────────────────────────────────

% Convert rotation degree to radians
rx_rad = deg2rad(objectRotationX_degree);
ry_rad = deg2rad(objectRotationY_degree);
rz_rad = deg2rad(objectRotationZ_degree);

% Rotation matrices
Rx_matrix = [1, 0,          0,         0;
             0, cos(rx_rad),  sin(rx_rad), 0;
             0, -sin(rx_rad), cos(rx_rad), 0;
             0, 0,          0,         1];

Ry_matrix = [cos(ry_rad), 0, -sin(ry_rad), 0;
             0,         1, 0,          0;
             sin(ry_rad), 0, cos(ry_rad),  0;
             0,         0, 0,          1];

Rz_matrix = [cos(rz_rad), sin(rz_rad), 0, 0;
        -sin(rz_rad), cos(rz_rad), 0, 0;
        0,          0,         1, 0;
        0,          0,         0, 1];

R_matrix = Rx_matrix * Ry_matrix * Rz_matrix;

% Scaling matrices
S_matrix = [Sx, 0,  0,  0;
            0,  Sy, 0,  0;
            0,  0,  Sz, 0;
            0,  0,  0,  1];

% Translation matrices       
T_matrix = [1, 0, 0, 0;
            0, 1, 0, 0;
            0, 0, 1, 0;
            Tx, Ty, Tz, 1];

% Calculate the World matrix
World = R_matrix * S_matrix * T_matrix;

% CAMERA INITIALIZATION AND CALCULATIONS
%─────────────────────────────────────────────────────────────────────────────────────────────

% Calculate where the camera is looking from
camera = [cameraPositionX, cameraPositionY, cameraPositionZ];
target = camLookAt;
worldUp = [0, 1, 0];

% Compute camera forward
forward = target - camera;
forward = forward / norm(forward);

% Compute camera right
right = cross(forward, worldUp);
right = right / norm(right);

% Compute camera up
up = cross(right, forward);

% rotation_look matrix for camera view
rotation_look = [ right(1),    up(1),    forward(1),   0;
                  right(2),    up(2),    forward(2),   0;
                  right(3),    up(3),    forward(3),   0;
                  0,           0,        0,        1 ];

% translation_look matrix for camera view
translation_look = [ 1, 0, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, 0;
                     -camera(1), -camera(2), -camera(3), 1 ];

% View matrix for camera
view = translation_look * rotation_look;

% PERSPECTIVE MATRIX INITIALIZATION AND CALCULATIONS
%─────────────────────────────────────────────────────────────────────────────────────────────

% Convert the field-of-view from degrees into radians and compute cotangent
a = deg2rad(fov);
cotHalf = 1 / tan(a/2);

% Perspective matrix
perspective_matrix = [ (1/aspect)*cotHalf,   0,                  0,               0;
                        0,                 cotHalf,              0,               0;
                        0,                   0,       farZ/(farZ - nearZ),        1;
                        0,                   0,  - (farZ * nearZ)/(farZ - nearZ), 0 ];

% CHANGING THE TRIANGLES FROM OBJECT-SPACE TO WORLD-SPACE
%─────────────────────────────────────────────────────────────────────────────────────────────

% Read the triangle data from the specified file
rawData = dlmread(filename);
[numTri, c] = size(rawData);

% Initialize an array to get the coordinates
verts_hom = zeros(numTri*3, 4);

% Iterate through the rows of rawData and create the triangle vertices
for i = 1:numTri
    row9 = rawData(i, :);
    v1 = row9(1:3);           
    v2 = row9(4:6);            
    v3 = row9(7:9);           

    % Convert to coordinates (x,y,z,1)
    baseIdx = (i-1)*3;              
    verts_hom(baseIdx+1, 1:3) = v1;
    verts_hom(baseIdx+2, 1:3) = v2;
    verts_hom(baseIdx+3, 1:3) = v3;
    verts_hom(baseIdx+(1:3), 4) = 1;
end

% Multiply the vertices by the World matrix
verts_world_hom = verts_hom * World;
verts_world = verts_world_hom(:, 1:3);

% USE BACKFACE CULLING TO REMOVE TRIANGLES FACING AWAY FROM THE CAMERA
%─────────────────────────────────────────────────────────────────────────────

% Iterate through each triangle to determine if it is facing the camera
keepMask_BF = false(numTri, 1);
for i = 1:numTri
    baseIdx = (i-1)*3;
    W1 = verts_world(baseIdx + 1, :);
    W2 = verts_world(baseIdx + 2, :);
    W3 = verts_world(baseIdx + 3, :);

    % Compute the edges and normal vector for backface culling
    e1 = W2 - W1;
    e2 = W3 - W1;
    N  = [e1(2)*e2(3) - e1(3)*e2(2), e1(3)*e2(1) - e1(1)*e2(3), e1(1)*e2(2) - e1(2)*e2(1)];

    % If dp > 0, the triangle is facing the camera thus keep it
    toCamera = [cameraPositionX, cameraPositionY, cameraPositionZ] - W1;
    dp = dot(N, toCamera);
    keepMask_BF(i) = (dp > 0);
end

% Keep track of the number of triangles kept after backface culling
numBFKept = nnz(keepMask_BF);
fprintf('Backface culling: Kept %d of %d triangles.\n', numBFKept, numTri);

% Reassemble all triangles that passed backface culling
rawData_BFkept = rawData(keepMask_BF, :); 
verts_hom_BFkept = zeros(3*numBFKept, 4);
idxBF = find(keepMask_BF);
for k = 1:numBFKept
    triIdx = idxBF(k);
    baseOrig = (triIdx - 1)*3;
    newBase = (k - 1)*3;
    verts_hom_BFkept(newBase + 1, :) = verts_hom(baseOrig + 1, :);
    verts_hom_BFkept(newBase + 2, :) = verts_hom(baseOrig + 2, :);
    verts_hom_BFkept(newBase + 3, :) = verts_hom(baseOrig + 3, :);
end

% Transform the culled mesh to camera space
verts_world_hom_BFkept = verts_hom_BFkept * World;   
verts_cam_hom_BFkept = verts_world_hom_BFkept * view;
verts_cam_BFkept = verts_cam_hom_BFkept(:, 1:3);

% Z-CLIPPING AND PERSPECTIVE PROJECTION
%─────────────────────────────────────────────────────────────────────────────

% Put camera-space vertices into clip space using the perspective matrix
verts_proj_hom = verts_cam_hom_BFkept * perspective_matrix;
w_coords = verts_proj_hom(:, 4);

% Convert the cordinates back by dividing by w
verts_ndc_all = [verts_proj_hom(:,1) ./ w_coords, verts_proj_hom(:,2) ./ w_coords, verts_proj_hom(:,3) ./ w_coords];

% Put the homogeneous coordinates back into the range of [-1, 1]
triPersp_all = zeros(numBFKept, 9);
for k = 1:numBFKept
    baseIdx = (k - 1)*3;
    v1 = verts_ndc_all(baseIdx + 1, :);
    v2 = verts_ndc_all(baseIdx + 2, :);
    v3 = verts_ndc_all(baseIdx + 3, :);
    triPersp_all(k,:) = [v1, v2, v3];
end

% Keep triangles within the viewing frustum else cull them out
keepMask_FR = false(numBFKept, 1);
for k = 1:numBFKept
    v1 = triPersp_all(k, 1:3);
    v2 = triPersp_all(k, 4:6);
    v3 = triPersp_all(k, 7:9);
    in1 = (abs(v1(1)) <= 1) && (abs(v1(2)) <= 1) && (v1(3) >= 0) && (v1(3) <= 1);
    in2 = (abs(v2(1)) <= 1) && (abs(v2(2)) <= 1) && (v2(3) >= 0) && (v2(3) <= 1);
    in3 = (abs(v3(1)) <= 1) && (abs(v3(2)) <= 1) && (v3(3) >= 0) && (v3(3) <= 1);
    keepMask_FR(k) = (in1 || in2 || in3);
end

% Tell how many triangles are kept after frustum culling
triPersp = triPersp_all(keepMask_FR, :);
numTri_final = size(triPersp, 1);
fprintf('Frustum culling after backface: Kept %d of %d triangles.\n', numTri_final, numBFKept);

% Map indices from backface-kept to final-frustum-kept for rawData
indicesBFkept = find(keepMask_BF);
indicesFR = find(keepMask_FR);
rawData_final = rawData_BFkept(keepMask_FR, :);

% Reassemble homogeneous vertices for the new set of triangles
verts_hom_final = zeros(3*numTri_final, 4);
for m = 1:numTri_final
    triIdx_BF = indicesBFkept(indicesFR(m));
    baseOrig  = (triIdx_BF - 1)*3;
    newBase   = (m - 1)*3;
    verts_hom_final(newBase + (1:3), :) = verts_hom(baseOrig + (1:3), :);
end

% Transform final triangles back into world and camera space for rendering
verts_world_hom_final = verts_hom_final * World;
verts_cam_hom_final = verts_world_hom_final * view;
verts_cam_final = verts_cam_hom_final(:, 1:3);

% Sort the the triangles based on which ones are closest to the camera and create from front to back
%─────────────────────────────────────────────────────────────────────────────

% Compute each triangle average Z value
avgZ = mean(triPersp(:, [3, 6, 9]), 2);
[~, sortIdx] = sort(avgZ, 'descend');
tri_sorted = triPersp(sortIdx, :);

% Build vertices_sorted and faces_sorted
vertices_sorted = reshape(tri_sorted.', 3, []).';
faces_sorted = reshape(1:(3*numTri_final), 3, []).';

% Reorder world-space vertex positions to match the sorted triangle order
verts_world_final = verts_world_hom_final(:, 1:3);
W_sorted = zeros(3*numTri_final, 3);
for newTriIdx = 1:numTri_final
    oldTriIdx   = sortIdx(newTriIdx);
    rowStart_old = (oldTriIdx - 1)*3 + 1;
    rowStart_new = (newTriIdx - 1)*3 + 1;
    W_sorted(rowStart_new:(rowStart_new + 2), :) = verts_world_final(rowStart_old:(rowStart_old + 2), :);
end

% Draw the final 3D → 2D graphics projection
%─────────────────────────────────────────────────────────────────────────────

% Set up the figure for rendering
figure('Name','3D → 2D Converter');
hold on;
axis equal;
xlim([-1 1]);
ylim([-1 1]);
xlabel('x\_axis');
ylabel('y\_axis');
title('3D → 2D Converter');

% Set up the triangles for the patch function
for t = 1:numTri_final
    % Get vertex indices for triangles
    faceIdx2D = faces_sorted(t, :);
    vertsThisTri2D = vertices_sorted(faceIdx2D, 1:2);  % (3×2)
    
    % Determine corresponding rows in W_sorted for world-space verts
    worldRows = (3*(t-1)+1) : (3*(t-1)+3);
    W1 = W_sorted(worldRows(1), :);
    W2 = W_sorted(worldRows(2), :);
    W3 = W_sorted(worldRows(3), :);

    % Compute world-space face normal
    e1 = W2 - W1;
    e2 = W3 - W1;
    Nw = cross(e1, e2);
    if norm(Nw) < 1e-8
        % If the light normal is too small, the triangle is dark
        faceColor = [0, 0, 0];
    else
        % Normalize normal and compute center
        n_hat   = Nw / norm(Nw);
        center = (W1 + W2 + W3) / 3;
        
        % Compute light direction vector and its length (strength)
        Lvec = lightPosition - center;
        L_norm = norm(Lvec);
        
        if L_norm < 1e-8
            lambert = 1;
        else
            % Normalize the light vector
            l_hat = Lvec / L_norm;
            lambert = max(0, dot(n_hat, l_hat));
        end
        % Determine the face color
        faceColor = M_diff * lambert;
    end

    % Render triangle in 2D with computed faceColor
    patch(vertsThisTri2D(:,1), vertsThisTri2D(:,2), faceColor, 'EdgeColor','none', 'FaceAlpha',1.0);
end