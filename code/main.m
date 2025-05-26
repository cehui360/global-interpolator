clear   
clc
tic
%% Input modeling_data
data = dlmread('modeling.txt');
minx = min(data(:,1));
maxx = max(data(:,1));
miny = min(data(:,2));
maxy = max(data(:,2));
%% Input checking_data
testdata = dlmread('checking.txt'); 
testdata(testdata(:,1)<minx,:) = [];
testdata(testdata(:,2)<miny,:) = [];
testdata(testdata(:,1)>maxx,:) = [];
testdata(testdata(:,2)>maxy,:) = [];
xt = testdata(:,1:2);
ze = testdata(:,3);
Num = size(testdata,1);
%% Parameter setting
para = struct('epsilon',0.001,'delta',10,'alpha',2,'beta',1,'max_iter',5);
% epsilon: parameter of the Γ-convergence  
% delta: smoothing parameter
%alpha/beta: regularization parameters
% max_iter: number of solution iterations
m = 1; % grid resolution
%% Grid point construction
xgrid = (minx - m/2):m:maxx;
ygrid = (miny - m/2):m:maxy;
[xg,yg] = meshgrid(xgrid,ygrid);
[r,c] = size(xg);
tablehead = struct('ncols',c,'nrows',r,'xllcorner',minx - m/2,'yllcorner',miny - m/2,'cellsize',m);
xg = reshape(xg,r*c,1);
yg = reshape(yg,r*c,1);
%% Compute and Define initial value
F = scatteredInterpolant(data(:,1:2),data(:,3),'natural','nearest');
zz1 = F(xg,yg); 
zs = reshape(zz1,r,c);
[r,c] = size(zs); 
s0 = ones(r,c);
z0 = ones(r,c);
u0 = zs;% initial values of u (g) 
%% Objective function solving 
[s, z, u] = VBCDA(s0, z0, u0, para,m,testdata,minx,miny,r,c);
u = reshape(u,r,c);
s = reshape(s,r,c);
z = reshape(z,r,c);
s = uint8(s);
z = uint8(z);
%u = uint8(u);
toc
zs = u;
%% Precision Validation
col = floor((testdata(:,1)-(minx-m/2))/m)+1;
row = floor((testdata(:,2)-(miny-m/2))/m)+1;
idxtest = (col-1)*r+row;
zz = zs(idxtest);
error0 = ze-zz;
error0(isnan(error0))=[];
rmse0 = sqrt(sum(error0.^2)/Num);
me0 = mean(abs(error0));
maxe0 = max(error0);
mine0 = min(error0);
fprintf(1,' MAE%.3fm RMSE%.3fm\n',me0,rmse0);
s1 = reshape(s,r*c,1);
z1 = reshape(z,r*c,1);
u1 = reshape(u,r*c,1);
u1 = [xg yg u1];  
idx1 = find(s1 == 0);
idx2 = find(z1 == 0);
break_points_s = u1(idx1,1:3);
break_points_z = u1(idx2,1:3);
%% Data Save
zs = flip(zs,1);
filename = strcat('griddem','.txt');
saveM2GisFile(filename,tablehead,zs)
dlmwrite('break_points_s', break_points_s, 'delimiter', '\t', 'precision', '%.6f');
dlmwrite('break_points_z', break_points_z, 'delimiter', '\t', 'precision', '%.6f');
s = flipud(s);
z = flipud(z);
%% Plot s\z
borderWidth = 1;
[rows, cols] = size(s);
s_with_border = zeros(rows + 2*borderWidth, cols + 2*borderWidth);
s_with_border(borderWidth+1:end-borderWidth, borderWidth+1:end-borderWidth) = s;
downsample_factor = 1;
s_with_border_small = s_with_border(1:downsample_factor:end, 1:downsample_factor:end);
figure('Color', [1 1 1]); 
imshow(s_with_border_small, 'InitialMagnification', 'fit', 'Colormap', [0 0 0; 1 1 1]);
axis equal;
axis off; 
title('s', 'Color', [0 0 0]); 
borderWidth = 1;
[rows, cols] = size(z);
z_with_border = zeros(rows + 2*borderWidth, cols + 2*borderWidth);
z_with_border(borderWidth+1:end-borderWidth, borderWidth+1:end-borderWidth) = z;
downsample_factor = 1;
z_with_border_small = z_with_border(1:downsample_factor:end, 1:downsample_factor:end);
figure('Color', [1 1 1]); 
imshow(z_with_border_small, 'InitialMagnification', 'fit', 'Colormap', [0 0 0; 1 1 1]);
axis equal; 
axis off; 
title('z', 'Color', [0 0 0]);
