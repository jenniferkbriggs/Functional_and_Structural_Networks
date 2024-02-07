%% Show network dot stick

%% This script shows the network on top of an islet image (e.g. Figure 6). 

%Jennifer Briggs 2020

close all 
clear all
clc


load('D:\SizeDependence\Analysis\Images\VidInfo.mat')
addpath('C:\Users\Jennifer Briggs\Documents\GitHub\Simulations');
addpath('C:\Users\Jennifer Briggs\Documents\GitHub\UniversalCode_Briggs\UniversalCode');
path = ['D:\SizeDependence\Analysis\Images\']
savename = ['D:\SizeDependence\Analysis\Images\']
Thr = .9
pic = 1;

folders = dir([path '*'])


for illy = 1:length(folders)
 Images = dir([path folders(illy).name '\*' '.lsm'])
if length(Images)>0
    
for lg = 6:length(Images)
   starttime = Vidinfo(illy).starttime(lg);
    endtime = Vidinfo(illy).endtime(lg);
%%{'39', '41','20','23', '1', '7', '10', '19','35', '43','22','24', '4', '11','14','15'}

number = lg

filepath = ['D:\SizeDependence\Analysis\Images\TrywFPS\' folders(illy).name '\' Images(lg).name]
imagepath = ['D:\SizeDependence\Analysis\Images\' folders(illy).name '\' Images(lg).name]

load([filepath '\CaWaveForm.mat'])

if pic

load([imagepath '\' 'Imaging.mat'])
filepath = ['D:\SizeDependence\Analysis\Images\TrywFPS\' folders(illy).name '\' Images(lg).name]
imagepath = ['D:\SizeDependence\Analysis\Images\' folders(illy).name '\' Images(lg).name]

   starttime = Vidinfo(illy).starttime(lg);
    endtime = Vidinfo(illy).endtime(lg);
% 
try
if zstacks == 1
    zz = 1
end
catch
    zstacks = 3
    zz = 2
    cachannel = 3
    howmanychannel = 3
end


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


if starttime == -1
    st=1;
else
    st = starttime;
end

if endtime == -1
    ed=length(T);
else
    ed=endtime;
end

T = T(st:ed);


images=double(IMG); % converts images to double precision
images = images(:,:,cachannel:howmanychannel:end);

RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable

tic
images = images(:,:,zz:zstacks:end);
images = images(:,:,st:ed-1);

sx=size(images,1);
sy=size(images,2);
sz=length(T);
for i=1:size(images,3)
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]); %applies filter to clean up images
    %nuimages(:,:,i)=medfilt2(nuimages(:,:,i),[5 5]);
end
toc
end


load([filepath '\' 'Masksv2.mat'])
load([imagepath '\' 'CellNumber.mat'])

% load([filepath '\' 'CellNumber.mat'])

if pic 
ImAv = sum(images,3); %compresses all frames into single array of intensities
HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
ImAvn = ImAv/max(ImAv(:));
HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
HSV(:,:,1) = 0.3333;%converts image to green image
RGB2 = hsv2rgb(HSV); %converts to rgb image
end
%%

for ll = 1:numcells
    [xp,yp] = find(CellMask == ll);
    x(ll) = mean(xp);
    y(ll) = mean(yp);
end

bad = find(isnan(x))

x(bad) = [];
y(bad) = [];
CellTC(:, bad) = [];
numcells = numcells - length(bad)
%%
fig = figure
if pic 
%RGB= im2gray(RGB2);
%RGB = imadjust(RGB);
imshow(RGB2)%, hold on
end

Adj = corr(CellTC);
Adj = Adj > Thr;
Adj = Adj - diag(diag(Adj));
AdjacencyGraph = graph(Adj);
Conn = sum(Adj);

Conarray = [0:max(Conn)];
Conarray2 = Conarray/length(Conarray)*100;
try
hubth = Conarray(Conarray2 > 60);
hubth = min(hubth);
Hubby = find(Conn >= hubth)
%Number of cells > 60% of Islet Links

    Nodec = repmat([.175 .54 .60],numcells,1);

for lll = 1:length(Hubby)
Nodec(Hubby(lll),:)=[.98 .122 .157];
end
catch
        Nodec = repmat([.175 .54 .60],numcells,1);

end
fig2 = figure
p = plot(AdjacencyGraph, 'Xdata',y,'YData',x, 'EdgeColor', 'b', 'NodeColor',Nodec,'MarkerSize',8, 'LineWidth',1 )
p.NodeLabel = [];
set(gca, 'YDir','reverse')
title([Images(lg).name(1:6) '  ' 'Mouse ' num2str(Vidinfo(illy).Mouse(lg)) ])
if Vidinfo(illy).type(lg) == 1
    check = contains(Images(lg).name(1:6), 'mm')
    M = 'mm'
elseif Vidinfo(illy).type(lg) == 2
    check = contains(Images(lg).name(1:6), 'pm')
    M = 'pm'
elseif Vidinfo(illy).type(lg) == 3
    check = contains(Images(lg).name(1:6), 'WT')
    M = 'pp'
end

if check ~= 1
    disp('Error')
end
set(gca, 'Visible', 'off')
set(gca, 'xticklabels', [])
set(gca, 'yticklabels', [])
% 
set(fig, 'Position', [100 100 1000 800])
set(fig2, 'Position', [100 100 1000 800])

try
saveas(fig, [path 'Images\FignNet\' M '\' 'BG' folders(illy).name '_' Images(lg).name '.png'])
saveas(fig2, [path 'Images\FignNet\' M '\' folders(illy).name '_' Images(lg).name '.png'])

catch
    mkdir([path 'Images\FignNet\' M ])
 saveas(fig, [path 'Images\FignNet\' M '\' 'BG' folders(illy).name '_' Images(lg).name '.png'])
saveas(fig2, [path 'Images\FignNet\' M '\' folders(illy).name '_' Images(lg).name '.png'])

end
% tit4 = title('D')
% tit4.FontSize = fs;
% ax = gca;
% ax.Position = [0.45 0 0.45 0.45];
clearvars x y
close(fig)
close(fig2)
end
end
end
