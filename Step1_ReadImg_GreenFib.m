
% This MATLAB SCRIPT reads a 3D confocal image stack of fibrin-platelet
% clot, where the first half of stack contains fibrin information in green,
% and second half of stack contains platelets in red. And converts the
% original image into extracted data with location, intensity of where
% red/green colors are detected.

% This script setup reads the fibrin in green. The output workspace of
% this script shoud be saved and then be used in step2-script.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
% cutoff = 190;
cutoff = 210;
dataIndex = 0;
filename = '2k-center_3.tif';
imageInfo = imfinfo(filename);
numberOfFrames = numel(imageInfo);
depth = numberOfFrames/2;

for z=1:depth
    
    originalImage = imread('2k-center_3.tif',z);
    originalImage = originalImage(400:end-220,200:end-220,:);
    dataIndex = dataIndex + 1
    image = originalImage(:,:,2);
%    image = adapthisteq(image,'NumTiles',[20, 20]);
     image = imsharpen(image,'Radius',50,'Amount',5);
%     
%     figure
%     imshow(image)
    cutoffImage = image;
       
    cutoffImage=imbinarize(cutoffImage,cutoff/255);
    cleanedImage = bwmorph(cutoffImage,'clean');
%     figure
%     imshow(cleanedImage)
%    cleanedImage = bwmorph(cleanedImage, 'bothat');
%    cleanedImage = bwmorph(cleanedImage, 'close');
    CC = bwconncomp(cleanedImage);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    indexes = (numPixels <= 100);

    toChange = {CC.PixelIdxList{indexes}};
    
        for a = toChange
            
            for b = 1:length(a{1})
                
                [i,j] = ind2sub(size(cleanedImage), a{1}(b));
                cleanedImage(i,j) = 0;
            
            end
            
        end
        
%    else
        
    
    %Average X and Y
    indexes = ~indexes;
    toCenters = {CC.PixelIdxList{indexes}};
    dilatedImage = bwmorph(cleanedImage,'dilate');
    var1 = numPixels(indexes)';
%     if dataIndex == 1
%         
%         numObjs = zeros(length(toCenters), depth);
%         numObjs(:,z) = numPixels(indexes)';
%         
%     end    
    
%   you can remove double for loop below
    for a = toCenters

        for b = 1:length(a{1})
            
            [i,j] = ind2sub(size(cleanedImage), a{1}(b));
            originalImage(i,j,:) = [255 255 255];

        end
%         scatter(averageX, averageY, 'o', 'filled')
%         pause(1/24);
    end

%    imshow(originalImage);
%     print(['layer' num2str(cutoff) '_z=' num2str(z) 'image'],'-dpng');
    
    sizes = size(dilatedImage);
    axis([0 sizes(2) 0 sizes(1)]);
    axis equal;
    
    hold on
%    pause(2)
    frameCount = 1;
    
    if dataIndex == 1
        
        numObjs = zeros(length(toCenters), depth);
        x = zeros(length(toCenters), depth);
        y = zeros(length(toCenters), depth);
        r = zeros(length(toCenters), depth);
        
    end
    
    index = 0;
    
    for a = toCenters
    
        index = index + 1;
        averageX = 0;
        averageY = 0;
        
        for b = 1:length(a{1})
        
            [i,j] = ind2sub(size(cleanedImage), a{1}(b));
            averageX = averageX + j/length(a{1});
            averageY = averageY + i/length(a{1});
        
        end
        
        numObjs(index, dataIndex) = var1(index);
        x(index, dataIndex) = averageX;
        y(index, dataIndex) = averageY;
        
        averageR = 0;
        
        for b = 1:length(a{1})
            
            [i,j] = ind2sub(size(cleanedImage), a{1}(b));
            Rx = j - averageX;
            Ry = i - averageY;
            R = ((Rx^2 + Ry^2)^(1/2));
            
            averageR = averageR + R/length(a{1});
        
        end
        
        r(index, dataIndex) = averageR;
            
    end

%     for i = 1:length(x)
%         
%         rad = r(i, dataIndex);
%         th = 0:pi/50:2*pi;
%         xunit = rad * cos(th) + x(i, dataIndex);
%         yunit = rad * sin(th) + y(i, dataIndex);
%         scatter(x(i, dataIndex), y(i, dataIndex), 'o', 'filled')
%         h = plot(xunit, yunit);
%               
%     end
%  
%        print(['objectCenters' num2str(cutoff) '_z=' num2str(z) 'image'],'-dpng');
%        frameCount = frameCount + 1;
%         
%        hold off

end

n = []; 
w = [];

for i = 1:(length(x(1,:))-1)
    
    if i ~= 1
        
        w{i-1} = n;
        n = [];
        
    end    
        
    for j = 1:length(x)
    
        for k = 1:length(x)
        
            if (x(j,i) - x(k,i+1)).^2 + (y(j,i) - y(k,i+1)).^2 < r(j,i).^2
            
                n = [n; i j k];
            
            end
            
        end
        
    end
    
end



w{i} = n;

k = 1;