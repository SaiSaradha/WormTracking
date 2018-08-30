function new_corners=SURF_BasedDetection(ref_img, image)

ref_img_gray = imsharpen(rgb2gray(ref_img));
%ref_img_gray=imfilter(ref_img_gray,h');

ref_pts = detectSURFFeatures(ref_img_gray);
[ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);

% Compare to video frame
%image = imread('/frames/136.jpg');
I = imsharpen(rgb2gray(image));
%I=imfilter(I,h');
% Detect features
I_pts = detectSURFFeatures(I);
[I_features, I_validPts] = extractFeatures(I, I_pts);
% figure;imshow(image);
% hold on; plot(I_pts.selectStrongest(50));

% Compare card image to video frame
index_pairs = matchFeatures(ref_features, I_features);

ref_matched_pts = ref_validPts(index_pairs(:,1)).Location;
I_matched_pts = I_validPts(index_pairs(:,2)).Location;

% figure, showMatchedFeatures(image, ref_img, I_matched_pts, ref_matched_pts, 'montage');
% title('Showing all matches');

% Define Geometric Transformation Objects
gte = vision.GeometricTransformEstimator; 
gte.Method = 'Random Sample Consensus (RANSAC)';

[tform_matrix, inlierIdx] = step(gte, ref_matched_pts, I_matched_pts);

ref_inlier_pts = ref_matched_pts(inlierIdx,:);
I_inlier_pts = I_matched_pts(inlierIdx,:);
% Transform the corner points 
% This will show where the object is located in the image

tform = maketform('affine',double(tform_matrix));
[width, height,~] = size(ref_img);
corners = [0,0;height,0;height,width;0,width];
new_corners = tformfwd(tform, corners(:,1),corners(:,2));
% figure;imshow(image);
% patch(new_corners(:,1),new_corners(:,2),[0 1 0],'FaceAlpha',0.5);
end
