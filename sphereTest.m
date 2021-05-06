clear
close all

load samplePoints pts

% Normal estimation
pcaNbrs = 20;

% Normal orientation
coarsegrid = 10.0; % in the units of the scan
graphNbrs = 5;

% Off-surface points
h = 1; % num off-surf points = 2 * numpts / h
L = 10.0; % height of off-surf pts

% Octree subdivision
Nmin = 2000;
Nmax = 5000;
expand = 1.1;

% RBF interpolation
rbf = @(r) r.^3;
polydegree = 2;
rho = 1e-5; % constant
regFn = @(N,t) 96*pi*N*t;

% Evaluation
alpha = 2*L;
Ngrid = 80; % for marching cubes
Hmax = L/2; % for tet mesh gen

% Turn waitbars on/off
global waitbarToggle
waitbarToggle = true;

%%

OBJ = surfFromPC(pointCloud(pts));

OBJ.estimateNormals(pcaNbrs);

OBJ.orientNormals(coarsegrid, graphNbrs);

OBJ.addOffSurfPoints(h, L);

OBJ.buildOctree(Nmin, Nmax, expand);

condition = OBJ.fit(rbf, polydegree, rho, regFn);

% For gridded sample & marching cubes:
[X, Y, Z, INTERP] = OBJ.gridSample(Ngrid, alpha);
S_cube = isosurface(X,Y,Z,INTERP);

% For tetrahedral mesh of alphaShape & marching tetra:
S_tet = aShapeSample(OBJ, alpha, Hmax);


%% Save

saveid = fix(clock);
saveid = sprintf('%d_%02d_%02d_%02d%02d/',saveid(1:5));

savepath = ['saves/', saveid];
mkdir(savepath)

save([savepath,'workspace.mat'])

OBJ.makeLog([savepath,'log.txt'])

%% Figure: Surface

fig = figure;
patch(S_cube, "LineStyle", "none", "FaceColor", "green")
view(3)
axis equal vis3d
camlight
title("Triangulation of the implicit surface by marching cubes")

savefig(fig, [savepath,'surface_cube'])

fig2 = figure;
patch(S_tet, "LineStyle", "none", "FaceColor", "green")
view(3)
axis equal vis3d
camlight
title("Triangulation of the implicit surface by marching tetrahedra")

savefig(fig2, [savepath,'surface_tet'])
