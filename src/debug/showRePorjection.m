function varargout = showRePorjection(pixels, blobs, imgSize)
%showRePorjection 粒子重投影和原图的关系
%   pixels 粒子重投影
img = zeros([imgSize, 3]);
im1 = zeros(imgSize); im1(pixels) = 1;
im2 = zeros(imgSize);
im3 = zeros(imgSize); 
for i = 1:1:length(blobs)
    idx1 = blobs(i).idx1;
    im2(idx1) = 1;
end
img(:,:,1) = im1;
img(:,:,2) = im2;
img(:,:,3) = im3;
if nargout == 0
	figure(1), imshow(img);
else
    varargout{1} = img;
end
end

