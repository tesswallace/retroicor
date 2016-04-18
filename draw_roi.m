function [ BW ] = draw_roi(image)
%Define region of interest
%returns [BW] a binary threshold of pixels to include in the stat image

h = figure;

imshow(image, []);    

title('Draw ROI...');

roi=impoly;

wait(roi);

BW=createMask(roi);
    
BW = double(BW);

close(h);

end

