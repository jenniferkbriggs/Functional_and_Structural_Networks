%% THIS PROGRAM IMPORTS .czi or .lsm CALCIUM IMAGING FILES FOR CELL-BY-CELL ANALYSIS OF CALCIUM TRACE TO IDENTIFY Cell Network
%% REAL-TIME USER INPUT REQUIRED.
%% Jennifer Briggs, Feb 2021
%% HOUSEKEEPING
close all
clear all
clc


cachannel = 1; %if there are multiple colors (e.g. calcium and nuclear stain, put them here)
nuchannel = 2;
howmanychannel = 1; %total number of channels (1 if only calcium)


filepath = %here you put the full path to the directory with images
ending = '.lsm';
savename = '/PlaywithSTD.mat'
numad = 1;
zstacks = 1;
zz = 1;
savepath = './ExampleData/'
ff = 1
starttime = -1
endtime = -1
gstart = 1


%here we find all of the imaging files and load them - this will need to be
%adjusted based on the user's computer and directory structure
for l = 1
    folderpaths = dir([filepath '*'])
    kn = []
    for illy = 1:length(folderpaths)
        if contains(folderpaths(illy).name, '.zip')
            kn = [kn illy];
        end
    end
    folderpaths(kn) = [];
    folderpaths([8 11]) = [];
    for g = gstart:length(folderpaths)
    starttime = Vidinfo(g).starttime;
    endtime = Vidinfo(g).endtime;
    filename = 'Analysis'
    F = dir([folderpaths(g).folder '/' folderpaths(g).name '/' folderpaths(g).name '/' '*' ending]);
    if length(starttime) > 0
    if length(F) ~= length(starttime)
        disp('files and file info do not match')
    end
    for ff = 1:length(starttime)
    savepath = ['./' folderpaths(g).name '/' F(ff).name];
    loadpath = ['.' folderpaths(g).name '/' F(ff).name];
    filename = [folderpaths(g).folder '/' folderpaths(g).name '/' folderpaths(g).name '/' F(ff).name]
    disp([folderpaths(g).name  '  ' F(ff).name])

    
%% LOADING THE CA IMAGE FILE
if exist('R') == 0
    try
        load([loadpath savename], 'R', 'zz', 'zstacks', 'cachannel', 'howmanychannel')
    catch
        R = bfopen([filename]); % Uses bfopen program to open .czi/.lsm image files
        try
            save([savepath savename],  '-V7.3')
        catch
            mkdir([savepath])
            save([savepath savename],  '-V7.3')
        end
     end
end


%% 
% omeMeta = 0
{1,4};
% voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
% voxelSizeXdouble = voxelSizeX.doubleValue();                                  % The numeric value represented by this object after conversion to type double                             
tic
pics=R{1};
pics=pics(:,1);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
pn = length(pics);

for i=1:pn
    IMG(:,:,i)=pics{i};
end

try
    for i=1:length(pics)
        T(i)=R{4}.getPlaneDeltaT(0, i-1).value;
    end
catch
    T=0:0.5:pn*0.5;
end
T = double(T);
T = T(cachannel:howmanychannel:end);


if starttime(ff) == -1
    st = 1;
else
    st = starttime(ff);
end

if endtime(ff) == -1
    ed=length(T);
else
    ed=endtime(ff);
end

T = T(st:ed);

images=double(IMG); % converts images to double precision
images = images(:,:,cachannel:howmanychannel:end);

RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable

clear pics R IMG;
% clear omeMeta;

output_dir = savepath;
toc

%% DECLARING IMAGE PROPERTIES
tic
images = images(:,:,st:ed-1);
sx=size(images,1);
sy=size(images,2);
sz=length(T)-1;
for i=1:sz
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]); %applies filter to clean up images
end

%Makes all frames 1fps
% fps =Vidinfo(g).fps(ff);
fps = 1;
 ct =1
for i = 1:fps:sz
    images(:,:,ct) = images(:,:,i);
    ct = ct+1;
end
images = images(:,:,1:ct-1);
%ImAv = sum(images,3); %compresses all frames into single array of intensities
ImAv = mean(images,3); 
HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
ImAvn = ImAv/max(ImAv(:));
HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
HSV(:,:,1) = 0.3333;%converts image to green image

RGB2 = hsv2rgb(HSV); %converts to rgb image

RGB= im2gray(RGB2);
RGB = imadjust(RGB);
OGFig = figure(1);
imshow(RGB)
toc
%% 

%%
%% MASKING ISLET FOR CELLS



%% User Draws ROIs around each cell within the islet
% ROIs are saved in "CellMask" array and called back to throughout analysis
CellMask = double(zeros(sx,sy));
numcells = 1;

% 
try %try to load in cell masks if already in loadpath
  
    load([loadpath '\Masks.mat'])
    
    load([loadpath '\CellNumber.mat'])
    
    ImAv = sum(images.*logical(CellMask),3); %compresses all frames into single array of intensities
    HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
    ImAvn = ImAv/max(ImAv(:));
    HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
    HSV(:,:,1) = 0.3333;%converts image to green image
    RGB2 = hsv2rgb(HSV); %converts to rgb image
    imshow(RGB2)   

catch
    NoSigFig = figure('Name','Draw ROIs Around Cells of Interest');
    imshow(RGB);
    keyboard
    k = 1;
    while k > 0
        disp('Draw ROIs Around Cells of Interest')
        ROIMask = imfreehand(); %User draws region around cell
        ROIMask = createMask(ROIMask); %Mask is created from drawn region
        CellMask = CellMask + ROIMask.*numcells; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
        CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
        UInp = input('Select Additonal ROI? 1=Yes, 0=No, 2=Redraw \n'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
        if UInp == 1 %User input to determine if another region is drawn
            k = k+1;
            numcells = numcells+1;
        elseif UInp == 2
            numcells = numcells;
            k = k;
        elseif UInp ~= 1
            k = 0;
        end
    end
    save([savepath '\Masks.mat'],'CellMask')
    save([savepath '\CellNumber.mat'],'numcells')
    close(NoSigFig);
    clear NoSigFig;  
end

%% CALCULATING/PLOTTING INTENSITY OVER TIME FOR EACH ROI
sz_2 = size(images,3);
PlotHandles = zeros(1,numcells);
PlotLabels = cell(1,numcells);
CellTC = zeros(sz_2,numcells);
TC = zeros(sz_2,1);

tic

%%Removes any pixels that do not belong to cell
CellTC = STDanalysis(images, CellMask);


% Plotting traces for entire time course
TCFig = figure('Name','Average Intensity Over Time');
plot(CellTC);
legend(PlotLabels);
title(F(ff).name)
%clear images MaskedIMGstack;
try
saveas(TCFig,[savepath '\Cellintestiy.tif']); %Saves figure of each cell's timecourse
save([savepath '\CaWaveForm.mat'],'CellTC')
catch
    mkdir(savepath)
    saveas(TCFig,[savepath '\Cellintestiy.tif']); %Saves figure of each cell's timecourse
end
save([savepath '\CaWaveForm.mat'],'CellTC')
numad = 1+ numad
clearvars -except Vidinfo starttime endtime zstacks howmanychannel nuchannel cachannel savename gstart ending filepath ff F imagepaths folderpaths numad AdjacencyAll l g zz zstacks cachannel
    end


        
