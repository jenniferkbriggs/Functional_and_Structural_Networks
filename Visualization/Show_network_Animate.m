%% This script shows the network as it changes over time. This scirpt was not used in the paper
close all 
clear all
clc

addpath('C:\Users\Jennifer Briggs\Documents\GitHub\Simulations');
addpath('C:\Users\Jennifer Briggs\Documents\GitHub\UniversalCode_Briggs\UniversalCode');
pic = 1;
    
filepath2 = ['D:\SizeDependence\Analysis\Images\10th Oct\Cx36mm42_55_11mMG2.lsm']
% F = dir([filepath2 '*.lsm'])
% 
% filepath2 = [filepath2 F(5).name]
if pic
load([filepath2 '\' 'CaWaveForm.mat'])

load([filepath2 '\' 'Imaging.mat'])


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
%%
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


images=double(IMG); % converts images to double precision
images = images(:,:,cachannel:howmanychannel:end);

RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable

tic
images = images(:,:,zz:zstacks:end);
st = 3000;
ed = 4000
T = T(st:ed);

sx=size(images,1);
sy=size(images,2);
sz=length(T);


for i=1:sz-1
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]); %applies filter to clean up images
    %nuimages(:,:,i)=medfilt2(nuimages(:,:,i),[5 5]);
end
toc
end
%% 
load([filepath2 '\' 'Masks.mat'])
load([filepath2 '\' 'CellNumber.mat'])

if pic 
ImAv2 = sum(images,3); %compresses all frames into single array of intensities
HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
mz = max(ImAv2(:));
ImAvn = ImAv2/max(ImAv2(:));
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
%%
figure

%%
for i = 21:sz-1

    ImAv = (images(:,:,i)); %compresses all frames into single array of intensities
HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
ImAvn = ImAv/max(ImAv(:));
HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
HSV(:,:,1) = 0.3333;%converts image to green image
RGB2 = hsv2rgb(HSV); %converts to rgb image
    if i < sz/3
% subplot('Position', [.02 .32 .45 .6])
% set(gca, 'xtick', [])
% set(gca, 'ytick', [])
set(gcf, 'PaperPosition', [0 0 4 2]);

% imshow(RGB2), hold on
% Adj = corr(CellTC(i-20:i,:));
% Adj = (Adj>.98);
% Adj = Adj-diag(diag(Adj));
% AdjacencyGraph = graph(Adj);
% plot(AdjacencyGraph, 'Xdata',y,'YData',x, 'EdgeColor', 'b', 'NodeColor','b','MarkerSize',.1, 'LineWidth',1 )
% %set(gca, 'YDir','reverse')
% title('Islet with Network')
% F1(i) = getframe(gcf) ;
% 
% drawnow
% 
% subplot('Position', [.51 .32 .45 .6])
set(gca, 'xtick', [])
set(gca, 'ytick', [])
imshow(RGB2)

% F1(i) = getframe(gcf) ;
drawnow
    if i == 21
        pause(5)
    end
else
imshow(RGB2), hold on
Adj = corr(CellTC(i-20:i,:));
Adj = (Adj>.96);
Adj = Adj-diag(diag(Adj));
AdjacencyGraph = graph(Adj);
plot(AdjacencyGraph, 'Xdata',y,'YData',x, 'EdgeColor', 'w', 'NodeColor','b','MarkerSize',.1, 'LineWidth',1,'NodeLabel', [])
%set(gca, 'YDir','reverse')
% 
% 
% subplot('Position', [.04 .1 .9 .18])
% plot(i,mean2(ImAv),'o')
% hold on
% xlim([21,sz])
% set(gca, 'xtick', [])
% set(gca, 'ytick', [])
% title('Average Intensity')
% F3(i) = getframe(gcf) ;
% 
drawnow

    end
end
 writerObj = VideoWriter([filepath2 'myVideo1.avi']);
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F1)
    % convert the image to a frame
    frame = F1(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


