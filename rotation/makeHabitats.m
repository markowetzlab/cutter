%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Clear previous variables
clear
close all

%%% Important variables to set manually:
cerr_path = '/Users/crispi01/CERR'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load CERR libraries
if isdir(cerr_path)
    addpath(genpath(cerr_path))
else
    error('Please specify the correct location of the CERR libraries.')
end

% Open image
cerrFileName = 'data/segmented_images.mat'
disp('Loading CERR file...')
planC = loadPlanC(cerrFileName,tempdir);
indexS = planC{end};

% Now reslice...
disp('Reslicing image...');
planC = reSliceScan(1,0.1,0.1,0.1,false,planC);
planC = reSliceScan(2,0.1,0.1,0.1,false,planC);
planC = reSliceScan(3,0.1,0.1,0.1,false,planC);
planC = reSliceScan(4,0.1,0.1,0.1,false,planC);
planC = reSliceScan(6,0.1,0.1,0.1,false,planC);
planC = reSliceScan(7,0.1,0.1,0.1,false,planC);
planC = reSliceScan(8,0.1,0.1,0.1,false,planC);

% Get masks
disp('Getting masks...')
[struc_3dprint] = getNonUniformStr(1, planC);
[struc_tum] = getNonUniformStr(3, planC);

% Create kidney
struc_rest = struc_3dprint;
struc_rest(struc_tum>0) = 0;

% Get the scan matrices
disp('Getting scan matrices...')
t1w_scan = getScanArray(planC{indexS.scan}(1));
d0_scan = getScanArray(planC{indexS.scan}(2));
perf_scan = getScanArray(planC{indexS.scan}(3));
hypercube_scan = getScanArray(planC{indexS.scan}(4));
ktrans_scan = getScanArray(planC{indexS.scan}(6));
t1map_scan = getScanArray(planC{indexS.scan}(7));
r2st_scan = getScanArray(planC{indexS.scan}(8));
scan_name = {'T1w', 'D0', 'PerfFrac', 'T2w', 'Ktrans', 'T1map', 'R2star'};

% Get scan arrays
z_array_scan_t1w = zscore(double(t1w_scan(struc_tum)));
z_array_scan_d0 = zscore(double(d0_scan(struc_tum)));
z_array_scan_perf = zscore(double(perf_scan(struc_tum)));
z_array_scan_hypercube = zscore(double(hypercube_scan(struc_tum)));
z_array_scan_ktrans = zscore(double(ktrans_scan(struc_tum)));
z_array_scan_t1map = zscore(double(t1map_scan(struc_tum)));
z_array_scan_r2st = zscore(double(r2st_scan(struc_tum)));
%
array_scan_t1w = double(t1w_scan(struc_tum));
array_scan_d0 = double(d0_scan(struc_tum));
array_scan_perf = double(perf_scan(struc_tum));
array_scan_hypercube = double(hypercube_scan(struc_tum));
array_scan_ktrans = double(ktrans_scan(struc_tum));
array_scan_t1map = double(t1map_scan(struc_tum));
array_scan_r2st = double(r2st_scan(struc_tum));

% Get x,y,z
indices = find(struc_tum);
[x,y,z] = ind2sub(size(struc_tum), indices);

% Create combined matrix
M = [x y z z_array_scan_t1w z_array_scan_d0 z_array_scan_perf z_array_scan_hypercube z_array_scan_ktrans z_array_scan_t1map z_array_scan_r2st];
M_notz = [x y z array_scan_t1w array_scan_d0 array_scan_perf array_scan_hypercube array_scan_ktrans array_scan_t1map array_scan_r2st];
Z = zscore(M);

% Apply clustering
rng(1)
disp('K-means clustering...')
num_clusters = [3];
col = brewermap(8,'Accent');

for k=num_clusters
    idx = kmeans(Z,k,'Replicates',20);

    figure;
    hold on;
    clusters_struct = zeros(size(struc_tum));
    clusters_struct(struc_tum) = idx;
    
    % Draw clusters
    clusNames = {}; 
    for eachClus=1:k
        new_struc = zeros(size(struc_tum));
        new_struc(struc_tum) = int8(idx==eachClus);
        blockPlot(new_struc, [0 0 0], 'FaceColor', col(eachClus,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        box on;
        grid on;
        titleString = sprintf('Cluster %i/%i', eachClus, k);
        title(titleString);
        clusNames{eachClus} = titleString;
    end
    
    % Draw kidney
    blockPlot(struc_rest, [0 0 0], 'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    clusNames{k+1} = 'Kidney';
    hold off;

    % Draw histograms
    figure;
    colours = {[220,245,166]./255, [239,129,245]./255, [249,157,155]./255};
    for pp=1:length(scan_name)
        subplot(3,3,pp);
        hold on;
        array_scan = M(:, pp+3);
        for ii=1:max(idx(:))
            histogram(array_scan(idx==ii),20,'FaceColor',colours{ii},'EdgeColor',[0.65,0.65,0.65]);
        end
        xlabel(scan_name{pp});
        box on;
    end
    legend({'Habitat 1','Habitat 2','Habitat 3'})
end

% Put the habitats back into the scan
scan_habitats = zeros(size(struc_tum));
scan_habitats(struc_tum) = idx;
scan_habitats(struc_rest) = num_clusters+1;
struc_merge = struc_tum | struc_rest;
bb_merge = getBox(struc_merge);
box_habitats = applyBox(scan_habitats,bb_merge);

% Apply the transformations...
load('transf_patient.mat')
rot1 = imwarp(box_habitats,tform1,'nearest');
rot2 = imwarp(rot1,tform2,'nearest');
rot3 = imwarp(rot2,tform3,'nearest');

% Get bounding vox of new volume
bb_rot3 = getBox(rot3);
box_final = applyBox(rot3,bb_rot3);

% Plot new, correctly-oriented volume
plotVolume(box_final,clusNames);

% Plot the slices
plotSlices(box_final,clusNames);



%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%

function inplanev = drawPlane(n,v,xv,yv,zv)
    inplanev = zeros(size(xv));
    offsets = [];
    for i=1:length(xv)
        x = xv(i);
        y = yv(i);
        z = yv(i);
        offset = abs( n(1)*(x-v(1)) + n(2)*(y-v(2)) + n(3)*(z-v(3)) );
        offsets(i) = offset;
        inplanev(i) = int16(offset<1);
    end
end

function [] = plotSlices(img, names)
    slices = 1:10:size(img,1);
    col = brewermap(8,'Accent');
    count = 1;
    for ii=1:length(slices)
        figure;
        theSlice = squeeze(img(slices(ii),:,:));
        rgb = label2rgb(theSlice,col);
        imshow(rgb,'InitialMagnification', 400)
        view([90 -90]);
        box on;
        grid on;
        hold on;
        title(['Slice ',num2str(count)]);
        dummyLegend(col, names, unique(theSlice));
        saveas(gcf,['slice',num2str(count,'%03.f'),'.png']);
        count = count+1;
        hold off
    end
end

function dummyLegend(c,names,elements)
    elements(elements==0) = []; % Remove zeros
    h = zeros(length(elements), 1);
    for k=1:length(elements)
        elem_k = elements(k);
        h(k) = plot(NaN,NaN,'o','MarkerFaceColor',c(elem_k,:),'MarkerEdgeColor',c(elem_k,:));
    end
    legend(h,names(elements));
end

function bb = getBox(mask)
ii = find(mask);
[u,v,w] = ind2sub(size(mask), ii);
u1 = min(u);
u2 = max(u);
v1 = min(v);
v2 = max(v);
w1 = min(w);
w2 = max(w);
bb = [u1,u2,v1,v2,w1,w2];
end

function box = applyBox(map,bb)
%// Now extract the bounding box
u1 = bb(1);
u2 = bb(2);
v1 = bb(3);
v2 = bb(4);
w1 = bb(5);
w2 = bb(6);
box = map(u1:u2, v1:v2, w1:w2);
end

function [] = plotVolume(thebox, names)
    figure;
    hold on;
    col = brewermap(8,'Accent');
    for i=1:length(names)
        s1 = zeros(size(thebox));
        s1(find(thebox==i)) = 1;
        blockPlot(s1, [0 0 0], 'FaceColor', col(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName',names{i});
    end
    legend
    box on
    grid on
    axis equal
    hold off;
end


function [fullMask] = getNonUniformStr(structNum, planC)
  scanNum = getStructureAssociatedScan(structNum, planC);
  indexS = planC{end};
  [rasterSegments, planC, isError] = getRasterSegments(structNum, planC);
  [mask3M, uniqueSlices] = rasterToMask(rasterSegments, scanNum, planC);
  scanArray3M = getScanArray(planC{indexS.scan}(scanNum));
  fullMask = false(size(scanArray3M));
  fullMask(:,:,uniqueSlices) = mask3M;
end


