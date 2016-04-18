function [ BW ] = draw_roi(image)
 
% Define region of interest
% This function requests the user to draw a region of interest using the 
% impoly tool and returns a binary mask [BW] of pixels to include in the 
% analysis

h = figure;

imshow(image, []);    

title('Draw ROI...');

roi=impoly;

wait(roi);

BW=createMask(roi);
    
BW = double(BW);

close(h);

end

