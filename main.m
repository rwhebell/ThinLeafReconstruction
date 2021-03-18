clear
close all

load samplePoints pts

% Denoising
denoiseNbrs = 50;
denoiseThresh = 1.0;

% Downsampling
method = 'gridaverage'; % choose gridaverage or random
downsampleParam = 3;
% for gridaverage: in the units of the scan (e.g., mm)
% for random: proportion between 0 and 1

% Normal estimation
pcaNbrs = 50;

% Normal orientation
coarsegrid = 3*downsampleParam; % in the units of the scan
graphNbrs = 10;

% Off-surface points
h = 1; % num off-surf points = 2 * numpts / h
L = 2*downsampleParam; % height of off-surf pts

% Octree subdivision
Nmin = 2000;
Nmax = 5000;
expand = 1.1;

% RBF interpolation
rbf = @(r) r.^3;
polydegree = 2;
% rho = [1e-10,10]; % range for GCV
rho = 1e-5; % constant
regFn = @(N,t) 96*pi*N*t;

% Evaluation
alpha = 2*L;
Ngrid = 100; % for marching cubes
Hmax = L/4; % for tet mesh gen

% Turn waitbars on/off
global waitbarToggle
waitbarToggle = true;

%%

OBJ = surfFromPC(pointCloud(pts));

OBJ.denoise(denoiseNbrs, denoiseThresh);

OBJ.downsample(method, downsampleParam);

OBJ.estimateNormals(pcaNbrs);

OBJ.orientNormals(coarsegrid, graphNbrs);

OBJ.addOffSurfPoints(h, L);

OBJ.buildOctree(Nmin, Nmax, expand);

condition = OBJ.fit(rbf, polydegree, rho, regFn);

% For gridded sample & marching cubes:
[X, Y, Z, INTERP] = OBJ.gridSample(Ngrid, alpha);
S = isosurface(X,Y,Z,INTERP);

% For tetrahedral mesh of alphaShape & marching tetra:
% S = aShapeSample(OBJ, alpha, Hmax);


%% Save

saveid = fix(clock);
saveid = sprintf('%d_%02d_%02d_%02d%02d/',saveid(1:5));

savepath = saveid;
mkdir(savepath)

save([savepath,'workspace.mat'])

OBJ.makeLog([savepath,'log.txt'])

%% Figure: Surface

fig = figure;
patch(S, "LineStyle", "none", "FaceColor", "green")
view(3)
axis equal vis3d
camlight

savefig(fig,[savepath,'surface'])

%% Figure: Cond number and smoothing parameter vs num points in patch

if length(rho)>1
	fig = figure; hold on
	
	yyaxis left
	scatter(OBJ.Ni, log10(condition), 16, 'LineWidth', 1)
	ylabel('log_{10}(condition number)')
	
	yyaxis right
	scatter(OBJ.Ni, log10(OBJ.fitOpts.rho), 30, 'sq', 'LineWidth', 1)
	ylabel('log_{10}(rho)')
	
	savefig(fig,[savepath,'condAndRho_vs_Ni'])
	close(fig)
end

%% Figure: Condition number vs smoothing parm

if length(rho)>1
	fig = figure;
	
	scatter(log10(OBJ.fitOpts.rho), log10(condition), ...
		[], OBJ.Ni, 'filled', 'MarkerEdgeColor', 'k' )
	xlabel('log_{10}(rho)')
	ylabel('log_{10}(condition number)')
	cb = colorbar;
	cb.Label.String = 'number in subdomain, N(i)';
	
	savefig(fig,[savepath,'cond_vs_rho'])
	close(fig)
end