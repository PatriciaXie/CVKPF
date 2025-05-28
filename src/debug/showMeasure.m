function showMeasure(img, blobs)
    im = zeros([size(img), 3]);
    im(:,:,2) = img;
    im(:,:,3) = img;
    for i = 1:1:length(blobs)
        idx = blobs(i).idx1;
        img(idx) = 1;
    end
    im(:,:,1) = img;
    figure, imshow(im);
end