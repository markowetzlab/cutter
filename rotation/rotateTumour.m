function [] = rotateTumour(cerr_path)
    % Load CERR libraries
    if isdir(cerr_path)
        addpath(genpath(cerr_path))
    else
        error('Please specify the correct location of the CERR libraries.')
    end

    % Open image
    cerrFileName = 'example_data/segmented_images.mat'
    disp('Loading CERR file...')
    planC = loadPlanC(cerrFileName,tempdir);
    indexS = planC{end};

    % Reslice image
    disp('Reslicing image...');
    planC = reSliceScan(1,0.1,0.1,0.1,false,planC);

    % Get masks
    disp('Getting masks...')
    [struc_inclfat] = getNonUniformStr(1, planC); % Everything, including fat
    [struc_tumour] = getNonUniformStr(3, planC);
    [struc_whole] = getNonUniformStr(4, planC); % Does not include fat
    [struc_edge1] = getNonUniformStr(5, planC);
    [struc_edge2] = getNonUniformStr(6, planC);
    [struc_contact1] = getNonUniformStr(7, planC);
    [struc_contact2] = getNonUniformStr(8, planC);
    [struc_contact3] = getNonUniformStr(9, planC);
    [struc_vessels] = getNonUniformStr(10, planC);
    [struc_hilum] = getNonUniformStr(11, planC);

    % Combine masks
    merged_mask = int16(struc_inclfat | struc_vessels | struc_edge1 | struc_edge2 | struc_contact1 | struc_contact2 | struc_contact3 );
    struc_combi = zeros(size(struc_inclfat));
    struc_combi(struc_inclfat) = 1;
    struc_combi(struc_whole) = 2;
    struc_combi(struc_hilum) = 8;
    struc_combi(struc_vessels) = 3;
    struc_combi(struc_edge1) = 4;
    struc_combi(struc_edge2) = 5;
    struc_combi(struc_tumour) = 6;
    struc_combi(struc_contact1) = 7;
    struc_combi(struc_contact2) = 7;
    struc_combi(struc_contact3) = 7;

    % Get bounding box
    bb_merged = getBox(merged_mask);
    box_combi = applyBox(struc_combi,bb_merged);

    % Plot original volume
    names = {'Fat','Whole','Vessels','Axis1','Axis2','Tumour','Contact','Hilum'};
    plotVolume(box_combi,names)

    % Apply 3 rotations
    [rot1,tform1] = rotateMidpoint(box_combi, 3, 7, 1, [1 2], [1 0 0]);
    [rot2,tform2] = rotateMidpoint(rot1, 3, 7, 1, [1 3], [0 0 1]);
    [rot3,tform3] = rotateCentroids(rot2, 4, 5, [1 2], [1 0 0]);

    % Get bounding vox of new volume
    bb_rot3 = getBox(rot3);
    box_final = applyBox(rot3,bb_rot3);

    % Plot new, correctly-oriented volume
    plotVolume(box_final,names);

    % Plot the slices
    plotSlices(box_final,names,'');

    % Save the transforms
    save('transf_patient.mat', 'tform1', 'tform2', 'tform3');

    % Smooth out and plot
    binary_box = zeros(size(box_final));
    binary_box(box_final~=0) = 1;
    smooth_box_final = smooth3(binary_box, 'gaussian', [9 9 9], 3);
    threshold_smooth_box_final = zeros(size(smooth_box_final));
    threshold_smooth_box_final(smooth_box_final>0.5) = 1;
    plotVolume(threshold_smooth_box_final,names);

    % Save the output file
    save('aligned_tumour.mat', 'box_final')
end


%%%%%%%%% AUXILIARY FUNCTIONS %%%%%%%%%%%

function [] = plotSlices(img, names, prefix)
    slices = 1:10:size(img,1);
    count = 1;
    num_parts = max(img(:));
    for ii=1:length(slices)
        figure;
        theSlice = squeeze(img(slices(ii),:,:));
        rgb = label2rgb(theSlice,jet(num_parts));
        imshow(rgb,'InitialMagnification', 500)
        view([90 -90]);
        box on;
        grid on;
        hold on;
        title(['Slice ',num2str(count)])
        dummyLegend(jet(num_parts), names, unique(theSlice))
        saveas(gcf,[prefix,'slice',num2str(count,'%03.f'),'.png'])
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
    %colors = {'r','g','b','k','y'}
    colors = jet(max(thebox(:)));
    for i=1:max(thebox(:))
        s1 = zeros(size(thebox));
        if (~any(thebox==i))
            continue
        end
        s1(find(thebox==i)) = 1;
        blockPlot(s1,[0 0 0], 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'DisplayName',names{i});
    end
    legend
    grid on
    box on
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


function [centroid] = getCentroid(img, num)
    mask = img==num;
    stats = regionprops(mask);
    centroid = stats.Centroid;
end

function [rotated1, tform] = rotateMidpoint(img, a_num, b_num, c_num, selD, x)
    % Get centroids
    a_centroid = getCentroid(img, a_num);
    b_centroid = getCentroid(img, b_num);
    c_centroid = getCentroid(img, c_num);

    % First rotation
    m_centroid = 0.5*(a_centroid+b_centroid);
    d = m_centroid-c_centroid;
    d1 = [0 0 0];
    d1(selD) = d(selD);
    d1 = d1./norm(d1);
    r = vrrotvec(d1,x);
    r = [r(1) r(2) r(3) pi-r(4)];
    m = vrrotvec2mat(r);
    m(:,4) = 0;
    m(4,:) = 0;
    m(4,4) = 1;
    tform = affine3d(m);
    rotated1 = imwarp(img,tform,'nearest');


    % Check rotation
    new_a_centroid = getCentroid(rotated1, a_num);
    new_b_centroid = getCentroid(rotated1, b_num);
    new_d = new_a_centroid-new_b_centroid;
    new_dnorm = new_d/norm(new_d);
end


function [rotated1, tform] = rotateCentroids(img, a_num, b_num, selD, x)
    % Get centroids
    a_centroid = getCentroid(img, a_num);
    b_centroid = getCentroid(img, b_num);

    % First rotation
    d = a_centroid-b_centroid;
    d1 = [0 0 0];
    d1(selD) = d(selD);
    d1 = d1./norm(d1);
    r = vrrotvec(d1,x);
    r = [r(1) r(2) r(3) pi-r(4)];
    m = vrrotvec2mat(r);
    m(:,4) = 0;
    m(4,:) = 0;
    m(4,4) = 1;
    tform = affine3d(m);
    rotated1 = imwarp(img,tform,'nearest');


    % Check rotation
    new_a_centroid = getCentroid(rotated1, a_num);
    new_b_centroid = getCentroid(rotated1, b_num);
    new_d = new_a_centroid-new_b_centroid;
    new_dnorm = new_d/norm(new_d);
end

function [binarybox] = addColumn(binarybox,filledbox,theNum)
    column_xs = false(size(filledbox,1), size(filledbox,2));
    % Map all voxels containing the structure onto the x-y plane
    for iz=1:size(filledbox,3)
        sel_voxels = squeeze(filledbox(:,:,iz)==theNum);
        column_xs(sel_voxels) = true;
    end

    % Smooth out cross section
    int_column_xs = zeros(size(column_xs));
    int_column_xs(column_xs) = 1;
    smooth_column_xs = imgaussfilt(int_column_xs, 3);
    bin_smooth_column_xs = false(size(column_xs));
    bin_smooth_column_xs(smooth_column_xs>0.5) = true;

    % Create column
    for iz=1:size(filledbox,3)
        this_slice = binarybox(:,:,iz);
        prev_values = this_slice(bin_smooth_column_xs);
        column_only_xs = bin_smooth_column_xs;
        column_only_xs(column_only_xs) = prev_values==0;
        this_slice(column_only_xs) = theNum;
        binarybox(:,:,iz) = this_slice;
    end
end
