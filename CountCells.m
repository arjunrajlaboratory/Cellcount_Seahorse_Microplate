function numNuclei=CountCells(imFilepath, plotOn,varargin)
ip=inputParser;
ip.addParameter('rescaleMinMax',[0 1000],@(x) isnumeric(x) && isequal(size(x),[1 2]) && (x(1)<x(2)))
ip.addParameter('junkSensitivity',0.4,@(x) isnumeric(x) && isscalar(x))
ip.addParameter('adaptThreshSensitivity',0.4,@(x) isnumeric(x) && isscalar(x))
ip.addParameter('minNucleusArea',15,@(x) isnumeric(x) && isscalar(x))
ip.addParameter('maxNucleusArea',200,@(x) isnumeric(x) && isscalar(x))
ip.addParameter('minNucleusCircularity',0.4,@(x) isnumeric(x) && isscalar(x))
ip.addParameter('gaussianBlurSigma',2,@(x) isnumeric(x) && isscalar(x))
ip.addParameter('junkDilatePx',15,@(x) isnumeric(x) && isscalar(x))

ip.parse(varargin{:}); % parse inputs from variable varargin
rescaleMinMax=ip.Results.rescaleMinMax;
junkSensitivity=ip.Results.junkSensitivity;
adaptThreshSensitivity=ip.Results.adaptThreshSensitivity;
minNucleusArea=ip.Results.minNucleusArea;
maxNucleusArea=ip.Results.maxNucleusArea;
minNucleusCircularity=ip.Results.minNucleusCircularity;
gaussianBlurSigma=ip.Results.gaussianBlurSigma;
junkDilatePx=ip.Results.junkDilatePx;

%%  read in image file
im=imread(imFilepath); % image is a uint16

% im=im(500:550,500:550); % for smaller test area
% show original image (showing images is optional)
if plotOn
    f=figure; clf;
    t=tiledlayout(2,6); t.TileSpacing='tight';
    nexttile
    imshow(im,[rescaleMinMax(1)*0.8 rescaleMinMax(2)*1.5]) % 0=black, 2000=white, for display purposes
    title('original image')
end
%% Rescale image
% rescale image (this will also convert it from uint16 to double, with values from 0 to 1)
lowIntensity=rescaleMinMax(1); % out of 65535 for uint16
highIntensity=rescaleMinMax(2);% out of 65535 for uint16

imRescaledDouble=rescale(im,'InputMin',lowIntensity,'InputMax',highIntensity); % to a double, 0 to 1

% optionally show this
if plotOn
    nexttile; linkaxes(t.Children)
    imshow(imRescaledDouble,[0 1])
    title(sprintf('imRescaled with low=%i, high=%i',lowIntensity,highIntensity))
end

imProc=imRescaledDouble;
%% Gaussian filter
useGaussianFilter=true;
if useGaussianFilter
    % Use gaussian filter to make image_blurred
    %gaussianBlurSigma=2;
    
    %%% First, we want to smooth the image and find regional max. These
    %%% will be the centroids of the cells.
    
    % set the gaussian filter to use:
    function_gaussian = @(block_struct) ...
        imgaussfilt(block_struct.data, gaussianBlurSigma);
    
    % smooth image with a gaussian (using block processing):
    image_blurred = blockproc(scale(imProc), [11000 11000], function_gaussian, 'BorderSize', [200 200]);
    imProc=image_blurred;
    
    % optionally show image_blurred
    if plotOn
        nexttile; linkaxes(t.Children)
        imshow(image_blurred,[0 1])
        title(sprintf('blurred with sigma=%i',gaussianBlurSigma))
    end
end
%% filter out big junk
junkFilterOn=true;
junkArea=10000;

if junkFilterOn
    %junkSensitivity=0.4;
    CCtemp=bwconncomp(imbinarize(imProc,adaptthresh(imProc, junkSensitivity)));
    rpTemp = regionprops(CCtemp,{'Area'});
    areaTemp = [rpTemp.Area]; % vector of area values
    idxJunk=areaTemp>=junkArea;
    PixelIdxListJunk=vertcat(CCtemp.PixelIdxList{idxJunk});
    junkImg=zeros(size(im),'logical'); % initialize
    junkImg(PixelIdxListJunk)=1;
    
    % dilate junkImg regions by a bit
    %junkDilatePx=15;
    se = strel('line',junkDilatePx,junkDilatePx); % x,y num dilation pixels
    junkImg = imdilate(junkImg,se); % %
    
    % show junk
    if plotOn
        nexttile; linkaxes(t.Children)
        imshow(junkImg,[],'Colormap',[0 0 0; 0.902, 0.902,0.9804])
        title(sprintf('junkImg with\njunkArea=%i,junkSensitivity=%0.2f\njunkDilatePx=%i', junkArea,junkSensitivity,junkDilatePx))
    end
    %
    %     mask=mask & ~junkImg;
    %
    %     nexttile; linkaxes(t.Children)
    %     imshow(mask)
    %     title(sprintf('mask after junk removal'))
    %imProc=imProc .* ~junkImg;
end

%% make a binary mask
useAdaptiveThreshold=true;
%adaptThreshSensitivity=0.4; % used if useAdaptiveThreshold=true; should be 0 to 1. In dentist2 it's 0.1 by default

if useAdaptiveThreshold
    T = adaptthresh(imProc, adaptThreshSensitivity); % for using an adaptive threshold
    mask = imbinarize(imProc,T); % for using an adaptive threshold, like is done in dentist2
    % show the adaptive thresh
    if plotOn
        nexttile; linkaxes(t.Children)
        imshow(T);
        title(sprintf('adaptive threshold\nwith sensitivity=%.2f',adaptThreshSensitivity));
    end
else
    mask = imbinarize(imProc);
end

mask=mask .* ~junkImg; % Added 15-Dec-2021

% show the mask
if plotOn
    nexttile; linkaxes(t.Children)
    imshow(mask)
    title(sprintf('mask with useAdaptiveThreshold=%i', useAdaptiveThreshold))
end



%% Use watershed
useWatershed=true;

if useWatershed
    % distance transform
    D = bwdist(~mask);
    if plotOn
        % show distance trasnform
        nexttile; linkaxes(t.Children)
        imshow(D,[prctile(D,0,'all') prctile(D,97,'all')])
        title('Distance Transform')
    end
    % calculate watershed transform
    L = watershed(-D); % on complement
    
    % make watershed label image only within the mask (Ie. candidate nuclei)
    L(~mask) = 0;
    
    % show watershed
    if plotOn
        nexttile; linkaxes(t.Children)
        rgb = label2rgb(L,'jet',[.5 .5 .5]);
        imshow(rgb)
        title('Watershed Transform')
    end
    
    % find connected components (contiguous '1' areas) and find properties on these
    CC = bwconncomp(L);
    rp = regionprops(CC,{'Area','Centroid','Circularity'});
    
else
    % find connected components (contiguous '1' areas) and find properties on these
    CC = bwconncomp(mask);
    rp = regionprops(CC,{'Area','Centroid','Circularity'});
end
%% choose parameter thresholds to filter out what is a real nucleus or not
%minNucleusArea=15;
%maxNucleusArea=200;
%minNucleusCircularity=.4; % circularity of a perfect circle is 1.

% then, determine which connected components match this filter
area = [rp.Area]; % vector of area values
circularity=[rp.Circularity];
idx = (area >= minNucleusArea) & (area <= maxNucleusArea) & (circularity>=minNucleusCircularity); % logical array where 1 means it meets these criteria

numNuclei=sum(idx);
fprintf('counted %i nuclei (of %i connected components) with useAdaptiveThreshold=%i sensitivity=%0.2f, minNucleusArea=%i, maxNucleusArea=%i, and minCircularity=%i \n',numNuclei,length(idx),useAdaptiveThreshold,adaptThreshSensitivity,minNucleusArea,maxNucleusArea,minNucleusCircularity)
%% Optionally, show image with centroids that meet these criteria
% first keep only centroids that meet these criteria:
% make a new rp (region props) structure that only includes real nuclei per
% these criteria
rpNuclei = rp(idx);

% get the centroids
centroids = [rpNuclei.Centroid];
centroids = reshape(centroids,2,[])'; % make it a 2-column array [rowInImage, columnInImage], where each row is another centroid

% show image with centroid (blue) overlays
if plotOn
    nexttile; linkaxes(t.Children)
    imshow(imRescaledDouble);
    hold on;% this makes it so that when we plot on the same axes, the next command doesn't delete the previous thing there (which in this case is the image)
    plot(centroids(:,1),centroids(:,2),'.','MarkerEdgeColor','blue')
    title(sprintf('%i good centroids (blue) with \nminArea=%i, maxArea=%i\n minCircularity=%.2f',numNuclei,minNucleusArea,maxNucleusArea,minNucleusCircularity))
    
    %% show good + bad (red) centroids together
    rpDiscarded=rp(~idx);
    centroidsBad=[rpDiscarded.Centroid];
    centroidsBad = reshape(centroidsBad,2,[])';
    
    nexttile; linkaxes(t.Children)
    imshow(imRescaledDouble); hold on;
    plot(centroids(:,1),centroids(:,2),'.','MarkerEdgeColor','blue');
    plot(centroidsBad(:,1),centroidsBad(:,2),'.','MarkerEdgeColor','red','MarkerSize',1)
    title('good (blue) and bad (red) centroids')
    
    %% Optionally overlay the area of the bad connected components
    PixelIdxListBad=vertcat(CC.PixelIdxList{~idx});
    badImg=zeros(size(im),'logical'); % initialize
    badImg(PixelIdxListBad)=1;
    
    PixelIdxListGood=vertcat(CC.PixelIdxList{idx});
    goodImg=zeros(size(im),'logical'); % initialize
    goodImg(PixelIdxListGood)=1;
    
    goodBadImg=zeros(size(im),'double'); % initialize
    goodBadImg(PixelIdxListGood)=1;
    goodBadImg(PixelIdxListBad)=0.5; % nothing=0, bad = 0.5, good=1;
    
    % plot it
    nexttile; linkaxes(t.Children)
    %imshow(imRescaledDouble)
    %hold on;
    %imshow(badImg,'Colormap',[0 0 0; 1 0 0])
    imshow(goodBadImg,[0 1],'Colormap',[0 0 0; 1 0.64 0; 0 1 1])
    title(sprintf('good (cyan) and bad (orange)\nconnected component areas'))
    hold on
    plot(centroids(:,1),centroids(:,2),'.','MarkerEdgeColor','blue');
    plot(centroidsBad(:,1),centroidsBad(:,2),'.','MarkerEdgeColor','red','MarkerSize',1)
    
end

end
