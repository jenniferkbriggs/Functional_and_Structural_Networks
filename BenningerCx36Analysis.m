%% THIS PROGRAM IMPORTS CALCIUM IMAGING FILES FOR CELL-BY-CELL ANALYSIS OF CALCIUM TRACE TO IDENTIFY Cell Network
%% REAL-TIME USER INPUT REQUIRED
%% Jennifer Briggs, Feb 2021
%% HOUSEKEEPING
close all
%clear all
clc
addpath('C:\Users\Jennifer Briggs\Documents\GitHub\Simulations');
addpath('C:\Users\Jennifer Briggs\Documents\GitHub\UniversalCode_Briggs\UniversalCode');
cachannel = 1;
nuchannel = 2;
howmanychannel = 1;
% Network opts
Opts.figs = 1;
Opts.printasSubplot = 0;
imagepath = struct;
load('/Volumes/Briggs_10TB/SizeDependence/Analysis/Images/VidInfo.mat')
filepath = '/Volumes/Briggs_10TB/SizeDependence/Analysis/Images/';
ending = '.lsm';
savename = '/PlaywithSTD.mat'
numad = 1;
zstacks = 1;
zz = 1;
savepath = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/STDanalysis/'
ff = 1
starttime = -1
endtime = -1
gstart = 1
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
    savepath = ['Volumes/Briggs_10TB/SizeDependence/Analysis/Images/TrywFPS/' folderpaths(g).name '/' F(ff).name];
    loadpath = ['Volumes/Briggs_10TB/SizeDependence/Analysis/Images/' folderpaths(g).name '/' F(ff).name];
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

if exist('R') == 0
try
   
    R = bfopen([filename]); % Uses bfopen program to open .czi/.lsm image files
    try
        save([savepath savename],  '-V7.3')
    catch
        mkdir([savepath])
        save([savepath savename],  '-V7.3')
    end
 end
end

%savepath = ['D:\Barak\Network\' imagepaths(l).name{g}];

% try
% if zstacks == 1
%     zz = 1
% end
% catch
%     zstacks = 3
%     zz = 2
% end

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
T = T(1:zstacks:end);
%T = T(3:3:end);

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
images = images(:,:,zz:zstacks:end);
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
%% MASKING ISLET

% User draws ROI around islet to plot whole islet timecourse
% User then determines where the start and end points for first responder
% and wave origin analyses
% imshow(RGB)
% disp('Draw ROI around entire islet');
% ROIMask_Start = imfreehand(); % User draws ROI around islet for average intensity
% ROIMask_Start = createMask(ROIMask_Start);
% sz_FR=ed-st+1;
% 
% StartMaskStack = images.*ROIMask_Start;
% IsletTC = shiftdim(mean(mean(StartMaskStack)),2);
% IsletTCfig = figure('Name','Whole Islet Time Course');
% plot(IsletTC); %plots intensity timecourse of entire islet so user can identify first responder and wave origin ranges
% %close(NoSigFig);
% clear StartMaskStack;
% % nuimages = nuimages.*ROIMask_Start;
% % %Getting rid of Background using JD code
% % nuHSV= ones(sx,sy,3); %preallocates a 3 dimensional array
% % nuHSV(:,:,3)= (mean(nuimages,3)./max(sum(nuimages,3))); %preallocates a 3 dimensional array
% % nuHSV(:,:,1) = 0.3333;%converts image to red image
% % RGBnu = hsv2rgb(nuHSV); %converts to rgb image
% 
% %% Getting rid of Background using JD code
% GrayFig = figure('Name','Remove Background');
% ImGray = rgb2gray(RGB); %converts image to gray
% imagesc(ImGray) %displays in default colormap; easier to see cell borders
% disp('Select Background Signal')
% bkgrndC=createMask(imfreehand); %allows user to select background area to be removed
% close(GrayFig);
% BackC=ImGray.*bkgrndC; %multiplies mask by the image
% BackC(~logical(BackC))=nan; %takes inverse of masked image and sets it to nan
% thresh = nanmax(nanmax(BackC)); %sets threshold based on the max of values after nan's are removed
% %close(OGFig);
% [ImgNon, ~] = RemovingAreaswNoSignal(ImGray,thresh); %Uses function from JD to remove background based on threshold
% ImGray(ImgNon)=0; %sets the removed areas to 0
% imagesc(ImGray);
% clear  GrayFig;



%% User Draws ROIs around each cell within the islet
% ROIs are saved in "CellMask" array and called back to throughout analysis
CellMask = double(zeros(sx,sy));
numcells = 1;

% 
try
    try
    load([loadpath '\Masksv2.mat'])
    catch
    load([loadpath '\Masks.mat'])
    end
    load([loadpath '\CellNumber.mat'])
    
    ImAv = sum(images.*logical(CellMask),3); %compresses all frames into single array of intensities
    HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
    ImAvn = ImAv/max(ImAv(:));
    HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
    HSV(:,:,1) = 0.3333;%converts image to green image
    RGB2 = hsv2rgb(HSV); %converts to rgb image

    imshow(RGB2)
%     check = input('Go on? 1 for yes or 2 for no')
%     if check == 2
%         return
%     end
        
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

% %Removes any pixels that do not belong to cell
% CellMasksave = CellMask;
% for i = 1:numcells
%     cctt = 1;
%     badpix = [2 2]; 
%     disp(i)
%     while ~isempty(badpix)
%     TCMask = CellMask; %Pulls in CellMask array
%     TCMask(TCMask ~= i) = 0; %Gets rid of all masks besides current one
%     MaskedIMGstack = images.*logical(TCMask); %Applies current mask to all frames of timecourse
%     %get rid of for loop
%     
%     [rr, cc] = find(mean(MaskedIMGstack,3));
%     for ii = 1:size(MaskedIMGstack,3)
%         TCnoZero = MaskedIMGstack(:,:,ii); %Pulls in current frame
%         TCnoZero = TCnoZero(TCnoZero>0); %Accounts for any zeros from preallocation
%         TCcheck(ii,:) = TCnoZero;
%         TC(ii) = mean(TCnoZero); %Calculates mean intensity of current frame
%     end
%      TCchecksm = smoothdata(TCcheck); %Smooth data
%     if size(TCchecksm,2) >0
%     for c=1:size(TCchecksm,2)
%             [C1, lag1] = xcorr(TCchecksm(:,c), median(TCchecksm'));
%             Maxcor(c) =  max(C1(find(lag1 > -2 & lag1 < 2)));       
%     end
%     badpix = find(Maxcor<1.2e6)
%     if length(badpix) > .75*size(TCcheck,2)       
%             badpix = find(Maxcor<mean(Maxcor) - 2*std(Maxcor))
%             cctt = cctt +1;
%     end
%    if cctt < 4  
% 
% %     figure, plot(rr,cc,'o')
%     for bb= 1:length(badpix)
%     %hold on, plot(rr(badpix(bb)),cc(badpix(bb)),'ro')
% 
%     CellMask(rr(badpix(bb)),cc(badpix(bb))) = 0;
%     end
%    else
% %        disp('Bad Pixels are too many')
% %        disp(savepath)
% %        fig = figure, plot(TCchecksm)
% %        x = input('Keep deleting or move on? 1 for move on, 0 for keep deleting')
% %        close(fig)
% %         if x == 0
% %             cctt = 1;
% %         else
% %             cctt = 5;
%             badpix = [];
% %         end
%    end
%         clear TCcheck Maxcor Maxcor2
%         
%         %figure, plot(TCchecksm)
%      else
%        disp('No pixels')
%        badpix = [];
%        clear TCcheck
%     end
%     end
%     
%     
%     CellTC(:,i) = TC; %Updates CellTC array with mean intensities from each frame of each cell
%     PlotLabels{i} = ['Cell' num2str(i)]; %Updates labels with current cell number for legend
% end
% try
% save([savepath '\Masksv2.mat'],'CellMask');
% catch
%     mkdir(savepath)
%     save([savepath '\Masksv2.mat'],'CellMask');
% end
% toc

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
%     end
%     end
% 
% end
% 

        
%%%%%%%%-----END OF CODE-----%%%%%%%%
