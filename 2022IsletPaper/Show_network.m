%% Show network dot stick
close all 
clear all
clc

addpath('C:\Users\Jennifer Briggs\Documents\GitHub\Simulations');
addpath('C:\Users\Jennifer Briggs\Documents\GitHub\UniversalCode_Briggs\UniversalCode');
Images = {'1', '7', '8','10', '19','39', '41', '23', '24', '3','4', '11', '14', '15','35', '43', '20', '22'} 
pic = 0;
for lg = 1:length(Images)
    
%%{'39', '41','20','23', '1', '7', '10', '19','35', '43','22','24', '4', '11','14','15'}
isletpath = Images{lg}

number = lg
load('', 'Network')
filepath = ['D:\Barak\Network\' Images{lg}]

if pic

load([filepath '\' 'Imaging.mat'])
filepath = ['D:\Barak\Network\' isletpath]


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

if isletpath == 24
    disp('l')
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


images=double(IMG); % converts images to double precision
images = images(:,:,cachannel:howmanychannel:end);

RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable

tic
images = images(:,:,zz:zstacks:end);

sx=size(images,1);
sy=size(images,2);
sz=length(T);
for i=1:sz
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]); %applies filter to clean up images
    %nuimages(:,:,i)=medfilt2(nuimages(:,:,i),[5 5]);
end
toc
end
if Images{lg}=='49'
load([filepath '\' 'Masks.mat'])
load([filepath '\' 'CellNumber.mat'])
else
try
load([filepath '\' 'Masksv2.mat'])
load([filepath '\' 'CellNumberv2.mat'])
catch
load([filepath '\' 'Masks.mat'])
load([filepath '\' 'CellNumber.mat'])
end
if pic 
ImAv = sum(images,3); %compresses all frames into single array of intensities
HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
ImAvn = ImAv/max(ImAv(:));
HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
HSV(:,:,1) = 0.3333;%converts image to green image
RGB2 = hsv2rgb(HSV); %converts to rgb image
end
%%
end
for ll = 1:numcells
    [xp,yp] = find(CellMask == ll);
    x(ll) = mean(xp);
    y(ll) = mean(yp);
end
%%
figure
if pic 
RGB= im2gray(RGB2);
RGB = imadjust(RGB);
imshow(RGB), hold on
end
Adj = Network(number).Rij;
Adj = Adj > .95;
Adj = Adj - diag(diag(Adj));
AdjacencyGraph = graph(Adj);
[Hubby] = find(((Network(number).N>round(size(Network(number).N,1)*.15)))/round(size(Network(number).N,1)))%Number of cells linked to more than 50% of Islet
Nodec = repmat([.175 .54 .60],numcells,1);
for lll = 1:length(Hubby)
Nodec(Hubby(lll),:)=[.98 .122 .157];
end
plot(AdjacencyGraph, 'Xdata',y,'YData',x, 'EdgeColor', 'b', 'NodeColor',Nodec,'MarkerSize',8, 'LineWidth',1 )
set(gca, 'YDir','reverse')
if lg < 10
    title('Control')
else
    title('KO')
end
set(gcf, 'Position', [100 100 1000 800])
saveas(gcf, ['D:\Barak\Network\Analysis\15_' num2str(isletpath) '.png'])
% tit4 = title('D')
% tit4.FontSize = fs;
% ax = gca;
% ax.Position = [0.45 0 0.45 0.45];

clearvars -except Images i isletpath lg pic
end