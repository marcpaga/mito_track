% Determines staining intensities lines from center to edge, calculates
% threshold based on set parameters.

[filename, pathname] = uigetfile({'*.*'},'Select file ...');
fullFileName = fullfile(pathname, filename);

prompt = {'Polarization starts at frame number:', ...
    'Polarization starting angle (0-360, noon is 90 degrees):', ...
    'Percentage below which staining is considered negative:', ...
    'Number of bins used for averaging angles:'};
dlg_title = 'Experiment variables:';
num_lines = 1;
defaultans = {'1','0','60','72'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

k1 = str2num(answer{1,1});
k1 = (k1 .* 2) - 1;
k2 = str2num(answer{2,1});
k3 = 0.01 * str2num(answer{3,1});
k4 = str2num(answer{4,1});

% determine edge and centroid of embryo based on myosin stain
myo.raw = imread(fullFileName, k1);
myo.bin = im2bw(myo.raw, graythresh(myo.raw));
myo.bin = bwareaopen(myo.bin, 100);
myo.label = bwlabel(myo.bin, 8);
center = regionprops(myo.label,'centroid');
myo.perim = bwperim(myo.bin);
[yperim, xperim] = find(myo.perim);

% load mitochondria image
mito.raw = imread(fullFileName, (k1 + 1));
meanmitointensity = mean2(mito.raw(myo.bin));

% because acosd angle determination spans 1 to 180 degrees rather than 1 to
% 360 --> split object into upper half and lower half, add 180 to lower.

upper = find(yperim < center.Centroid(2));
lower = find(yperim > center.Centroid(2));

% determine mito < threshold intensity edge for upper half
for uu = 1 : length(upper);
    x(uu) = xperim(upper(uu));
    y(uu) = yperim(upper(uu));
        
    [linex, liney] = bresenham(center.Centroid(1), ... 
        center.Centroid(2), x(uu), y(uu));

    for i = 1 : size(linex, 1);
        lineint(i) = mito.raw(liney(i), linex(i));
    end

    th1 = k3 * meanmitointensity;
    k = lineint > th1;
    seq = ones(1, 10);
    th2 = findstr(k, seq);
        
    if isempty(th2)
       upperxmitothreshold(uu) = x(uu);
    else
       upperxmitothreshold(uu) = linex(th2(end)+9);
    end
    
    if isempty(th2)
       upperymitothreshold(uu) = y(uu);
    else
       upperymitothreshold(uu) = liney(th2(end)+9);
    end
    
    upperdistancetoedge(uu) = sqrt((x(uu) - upperxmitothreshold(uu)).^2 + ...
        (y(uu) - upperymitothreshold(uu)).^2);
    
    % also calculate vector angles compared to three o'clock
    deltax(uu) = x(uu) - center.Centroid(1);
    deltay(uu) = y(uu) - center.Centroid(2);
    vectors(uu) = sqrt((deltax(uu).^2) + (deltay(uu).^2));
    upperangles(uu) = acosd(deltax(uu)/vectors(uu));
        
    % store line data
    LineStats(uu).Angles = upperangles(uu);
    LineStats(uu).MitoIntensities = lineint;
    LineStats(uu).ThresholdCoordinates(1) = upperxmitothreshold(uu);
    LineStats(uu).ThresholdCoordinates(2) = upperymitothreshold(uu);
    LineStats(uu).TotalLength = vectors(uu);
    LineStats(uu).ThresholdLength = upperdistancetoedge(uu);
    
    % clear lineint as otherwise its size seems fixed
    clear lineint i
end

% repeat mito < threshold intensity edge for lower half, add 180 to angles
for ul = 1 : length(lower);
    x(ul) = xperim(lower(ul));
    y(ul) = yperim(lower(ul));
        
    [linex, liney] = bresenham(center.Centroid(1), ... 
        center.Centroid(2), x(ul), y(ul));

    for i = 1 : size(linex, 1);
        lineint(i) = mito.raw(liney(i), linex(i));
    end

    th1 = k3 * meanmitointensity;
    k = lineint > th1;
    seq = ones(1, 10);
    th2 = findstr(k, seq);
    
    if isempty(th2)
       lowerxmitothreshold(ul) = x(ul);
    else
       lowerxmitothreshold(ul) = linex(th2(end)+9);
    end
    
    if isempty(th2)
       lowerymitothreshold(ul) = y(ul);
    else
       lowerymitothreshold(ul) = liney(th2(end)+9);
    end
    
    lowerdistancetoedge(ul) = sqrt((x(ul) - lowerxmitothreshold(ul)).^2 + ...
        (y(ul) - lowerymitothreshold(ul)).^2);
    
    % also calculate vector angles compared to three o'clock
    deltax(ul) = x(ul) - center.Centroid(1);
    deltay(ul) = y(ul) - center.Centroid(2);
    vectors(ul) = sqrt((deltax(ul).^2) + (deltay(ul).^2));
    lowerangles(ul) = acosd(deltax(ul)/vectors(ul));
    anglecorrection = 180 - lowerangles;
    correctedlowerangles = anglecorrection + 180;
    
    % store line data
    LineStats(uu+ul).Angles = correctedlowerangles(ul);
    LineStats(uu+ul).MitoIntensities = lineint;
    LineStats(uu+ul).ThresholdCoordinates(1) = lowerxmitothreshold(ul);
    LineStats(uu+ul).ThresholdCoordinates(2) = lowerymitothreshold(ul);
    LineStats(uu+ul).TotalLength = vectors(ul);
    LineStats(uu+ul).ThresholdLength = lowerdistancetoedge(ul);
    
    % clear lineint as otherwise its size seems fixed
    clear lineint i
end

% concatenate threshold coordinates
xthreshold = horzcat(upperxmitothreshold, lowerxmitothreshold);
ythreshold = horzcat(upperymitothreshold, lowerymitothreshold);
distancetoedge = horzcat(upperdistancetoedge, lowerdistancetoedge);
angles = horzcat(upperangles, correctedlowerangles);

% re-order angles based on k2: ...-1...k2...1...
angles = angles - k2;

for i = 1 : length(angles);
    if angles(i) < 0;
       angles(i) = angles(i) + 360;
    end
end

for j = 1 : length(angles);
    if angles(j) > 180;
       angles(j) = angles(j) - 360;
    end
end

% calculate distance to edge as percentage
for w = 1 : length(LineStats);
    LineStats(w).ThresholdPercentage = LineStats(w).ThresholdLength ./ ...
        LineStats(w).TotalLength * 100;
    thresholdpercentage(w) = LineStats(w).ThresholdLength ./ ...
        LineStats(w).TotalLength * 100;
    thresholdlength(w) = LineStats(w).ThresholdLength;
end

% create some check plots
imshow(myo.bin);
hold on
scatter(xperim, yperim);
hold on
scatter(center.Centroid(1), center.Centroid(2), 'filled', ... 
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
hold off

figure
imshow(mito.raw);
hold on
scatter(xperim, yperim);
hold on
scatter(xthreshold, ythreshold);
hold off

figure
scatter(angles, distancetoedge);
xlabel('line angle (degrees)')
ylabel('distance between edge and threshold (pixels)')

% average results over 24 bins
nbins = k4;
degs = 360 ./ nbins;
negdegs = -1 * degs;

for b1 = 1 : (nbins * 0.5);
    inbin = angles > (degs * (b1-1)) & angles < (degs * b1);
    bindistancetoedge = thresholdpercentage(inbin);
    binthresholdlength = thresholdlength(inbin);
    meandistancetoedge(b1, 1) = degs * (b1-1);
    meandistancetoedge(b1, 2) = degs * b1;
    meandistancetoedge(b1, 3) = ((degs * (b1-1)) + (degs * b1)) ./ 2;
    meandistancetoedge(b1, 4) = mean(bindistancetoedge);
    meandistancetoedge(b1, 5) = mean(binthresholdlength);
    BinStats(b1).BinRange(1) = meandistancetoedge(b1, 1);
    BinStats(b1).BinRange(2) = meandistancetoedge(b1, 2);
    BinStats(b1).ThresholdPercentage = meandistancetoedge(b1, 4);
    BinStats(b1).AverageLengthInPixels = meandistancetoedge(b1, 5);
end

for b2 = 1 : (nbins * 0.5);
    inbin = angles < (negdegs * (b2-1)) & angles > (negdegs * b2);
    bindistancetoedge = thresholdpercentage(inbin);
    binthresholdlength = thresholdlength(inbin);
    meandistancetoedge(b1+b2, 1) = negdegs * (b2-1);
    meandistancetoedge(b1+b2, 2) = negdegs * b2;
    meandistancetoedge(b1+b2, 3) = ((negdegs * (b2-1)) + (negdegs * b2)) ./ 2;
    meandistancetoedge(b1+b2, 4) = mean(bindistancetoedge);
    meandistancetoedge(b1+b2, 5) = mean(binthresholdlength);
    BinStats(b1+b2).BinRange(1) = meandistancetoedge(b1+b2, 1);
    BinStats(b1+b2).BinRange(2) = meandistancetoedge(b1+b2, 2);
    BinStats(b1+b2).ThresholdPercentage = meandistancetoedge(b1+b2, 4);
    BinStats(b1+b2).AverageLengthInPixels = meandistancetoedge(b1+b2, 5);
end

figure
scatter(meandistancetoedge(:,3), meandistancetoedge(:,4));
xlabel('bin average of line angles (degrees)')
ylabel('average distance between edge and threshold (percetage of line length)')

figure
scatter(meandistancetoedge(:,3), meandistancetoedge(:,5));
xlabel('bin average of line angles (degrees)')
ylabel('average distance between edge and threshold (pixels)')

% clear lose variables
clear deltax deltay i k lineint linex liney lower lowerangles vectors
clear correctedlowerangles lowerdistancetoedge seq th1 th2 ul upper 
clear upperangles upperdistancetoedge uu x y anglecorrection
clear lowerxmitothreshold upperxmitothreshold xthreshold ythreshold
clear lowerymitothreshold upperymitothreshold angles distancetoedge
clear nbins degs inbin w center b1 b2 bindistancetoedge binthresholdlength
clear xperim yperim k1 k2 prompt negdegs meanmitointensity num_lines
clear thresholdpercentage k3 k4 dlg_title answer defaultans j pathname
clear filename fullFileName meandistancetoedge thresholdlength