function [debut, fin, peak, start_node, end_node, scalingfactor] = aisWW3D(coordfile,imgfile,axonchan,thresh,varargin)

% September2009  assumes that AIS channel is last channel line 34
%%AISWW3D analyses the given image stack to determine where the AIS is.
%
% [DEBUT, FIN, PEAK] = AISWW3D(COORDFILE,IMGFILE,AXONCHAN,THRESH)
%
% INPUTS
% COORDFILE - full path to a text file with x,y,z positions of a manual 
% drawing along the axon. Use ImageJ's/Fiji's Simple Neurite Tracer to get 
% the trace.
% IMGFILE - full path to the stack you want to analyse, in tif format
% AXONCHAN - the number of the channel that contains the staining for the
% AIS
% THRESH - threshold above which the AIS staining is considered to be part
% of the actual AIS. Set to 0.25 for my analysis, Matt used 0.33.
%
% OUTPUTS
% DEBUT - start of the AIS in um
% FIN - end of the AIS in um
% PEAK - peak of the AIS in um
%
% See also AISWW3DMASTER, AISWW3DPLOT

    %%% open image file
    pic = imread(imgfile,1);              % read first image in the stack
    info = imfinfo(imgfile);              % all we know about the stack
    n = length(info);                     % number of images in stack
    description = info(1).ImageDescription;
    pixmic = info.XResolution;
    scalingfactor =  1/pixmic;
    
    [s e] = regexp(description,'channels=[0-9]'); %Alex - Index Location of Channels
    channelnum = str2double(description(s+9:e)); % number of channels
    axonchan = channelnum %% September2009  assumes that AIS channel is last channel
    
    [s e] = regexp(description,'spacing=[0-9].[0-9]');
    scalez = str2double(description(s+8:e)); % z-spacing
    if isnan(scalez) %If it doesnt find spacing infor assumes 1 um
        scalez = 1;
        warning('AISWW3D:NoZSpacing','No z-spacing information, using 1 um as default value.');
    end
    
    % NB second image is first image of second channel, third image is
    % first image of third channel, fourth image is second image of first
    % channel, etc.
    
    stack = double(zeros([size(pic),1,n]));  % empty 4D array to store the stack
    
    for i = 1:n                           % read the images into the array (Alex-All channels included)
        nextpic = double(imread(imgfile,i));
        stack(:,:,i) = nextpic;
    end
    
    %%% read 3D co-ordinates
    [x,y,z] = textread(coordfile,'%*d%*d%f%f%f%*[^\n]','headerlines',1);
    
    scalex = abs(x(1)-x(2));    % alternatively, use info(1).XResolution
    scaley = abs(y(1)-y(2));
    scale = max(scalex,scaley); %one of them might be zero, if only one coord changed in first step
    if scale ==0
        if exist('info(1).XResolution')
            scale = 1/info(1).XResolution;
        else
            scale = 0.12;
            warning('AISSYNAPSEDIL2:NoXYSpacing','No xy-spacing information, using 0.12 um as default value.');
        end
    end

    points = [y./scale+1,x./scale+1,z./scalez+1]; % add one to all points because ImageJ starts with 0, Matlab with 1
    points = round(points);  %for some reason, the z-coordinates all have some decimal places
    [~,I] = unique(points,'rows');
    points = points(sort(I),:); % put points back into original order
    
    if length(varargin) ==1
        points = varargin{1};
    end
    
    axonchanindex = axonchan:channelnum:n; %Alex-Index frames that belong to AIS channel
    axonchannel = stack(:,:,axonchanindex); % all frames that belong to the channel with the axon - Alex substack with only channel of interest (all rows, all columns, only z of interests)
    
    %%% detect AIS along the given path (i.e. axon), see ais.m line 334
    for i = 1:size(axonchannel,3) %1:Size of Z
        axonchannel(:,:,i) = nanmedfilt2(axonchannel(:,:,i));    %maybe smooth in z? only if z steps are small enough
    end
    
    axonvalues = axonchannel(sub2ind(size(axonchannel),points(:,1),points(:,2),points(:,3))); %Alex -selects just axon
    % now also smooth along axon, Matt uses sliding mean of 41 (20 each
    % side)
    for i = 1:length(axonvalues)
        d = 20; % how big a sliding window? if d=20, window is 41 pixels large 
        if i<(d+1)  % not enough pixels at the start
            axonsmooth(i) = mean([axonvalues(1:i); axonvalues(i+1:i+d)]); 
        elseif i>(length(axonvalues)-(d+1)) % not enough pixels at the end
            axonsmooth(i) = mean([axonvalues(i-d:i-1); axonvalues(i:length(axonvalues))]);
        else
            axonsmooth(i) = mean([axonvalues(i-d:i); axonvalues(i+1:i+d)]);
        end
    end
    
    %%% background drawing and calculation
%     axonimg = mean(axonchannel,3);
%     if isempty(bg)
%         roi_info = imSelectROI(mat2gray(axonimg));
%         bgrect.(strcat('img_',imgfile)) = [roi_info.Ymin, roi_info.Xmin, roi_info.DY, roi_info.DX];
%     end
%     rect = bgrect.(strcat('img_',imgfile));
%     bgarea = axonimg(rect(1):rect(1)+rect(3),rect(2):rect(2)+rect(4));
%     bgvalue = mean(bgarea(:));
        
    norm_axon = (axonsmooth - min(axonsmooth)) ./ (max(axonsmooth)-min(axonsmooth));
    %norm_axon = (axonsmooth - bgvalue) ./ (max(axonsmooth) - bgvalue);
    plot(norm_axon) % show the axon's intensity profile
    
    max_i = find(norm_axon==1);
    max_i = max_i(1);   % find the first peak/max value (starting from soma)
    pix_narray = (1:length(points));
    
    ais_end = find( (pix_narray>max_i) & (norm_axon<thresh));
    
    if ~isempty(ais_end)
        ais_end = ais_end(1); % point index along axon past max where fluorescence intensity falls to f of its peak,
    else                      % or just the last pixel along the line
        ais_end = length(points);
        warning('AISWW3D:AIS_too_short','End of AIS set to end of drawing. True end cannot be determined in: %s',coordfile)
    end
    ais_start = find( (pix_narray<max_i) & (norm_axon<thresh)); 
    if ~isempty(ais_start)
        ais_start = ais_start(length(ais_start));
    else
        ais_start = 1; % point index along axon pre max where fluorescence intensity falls to f of its peak
    end
    
    
 start_node = ais_start
 end_node = ais_end 
    
    %%% stretch AIS using the wisdom of Pythagoras of Samos
    
    stretchAIS = ones(size(points,1),1);
    for i = 1:size(points,1)-1
    
        p1 = points(i,:);
        p2 = points(i+1,:);
        dist = ([scale scale scalez].*(p2-p1)).^2;
        stretchAIS(i) = sqrt(sum(dist));
        
    end
    stretchAIS = cumsum(stretchAIS);
    
    if ais_start ==1
        debut = 0;
    else
        debut = stretchAIS(ais_start-1);  %%%% AIS start position in um
    end
       
 
    
%     [scs,p] = csaps(stretchAIS,norm_axon,0.9);
%     figure
%     fnplt(scs,2)
%     
%     [scs,p] = csaps(stretchAIS,norm_axon,0.85);
%     figure
%     fnplt(scs,2)
%     
%     [scs,p] = csaps(stretchAIS,norm_axon,0.8);
%     figure
%     fnplt(scs,2)
    
    [~,peak] = max(norm_axon(ais_start:ais_end));
    peak = stretchAIS(peak+ais_start-1);
    fin = stretchAIS(ais_end-1);      %%%% AIS end position in Ch1, in um
    lngth = fin-debut;   %%%% AIS length in Ch1, in um
    
    plot(stretchAIS,norm_axon)
    hold on
    plot(debut,0,'r*')
    plot(fin,0,'r*')
    title(imgfile)
    hold off
    

  
    
    