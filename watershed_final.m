img=imread('/frames/299.jpg');
I=imcrop(img);
I1=rgb2gray(I);
I2=im2bw(I1);
%I2=I1 <= 110;
D=bwdist(I2);
imshow(D)
L=watershed(-D);
imshow(L)
w=L==0;
figure, imshow(w);  
g2= I2 & ~w;
figure, imshow(g2);

%segmentation through boundaries:

bound_seg=bwboundaries(~I2,4,'holes');
figure, imshow(I2);
hold on;
for k = 1:length(bound_seg)
   b_curr = bound_seg{k};
   plot(b_curr(:,2), b_curr(:,1), 'LineWidth', 2); hold on;
end



%im1=frame_curr;
im2=frame_curr;
%bwim1=adaptivethreshold(im1,11,0.03,0);
bwim2=adaptivethreshold(im2,15,0.02,0);
%subplot(2,2,1);
%imshow(im1);
%subplot(2,2,2);
figure, imshow(bwim1);
%subplot(2,2,3);
%imshow(im2);
%subplot(2,2,4);
figure, imshow(bwim2);

%Adaptive thresholding for segmentation : Useful tool
figure, imshow(bwim2); hold on;
colors=['b' 'g' 'r' 'c' 'm' 'y'];
for k = 1:length(BBT)
boundary = BBT{k}; cidx = mod(k,length(colors))+1;
plot(boundary(:,2), boundary(:,1), colors(cidx), 'LineWidth', 2);
drawnow;
%randomize text position for better visibility
rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
col = boundary(rndRow,2); row = boundary(rndRow,1);
h = text(col+1, row-1, num2str(L(row,col)));
set(h,'Color',colors(cidx),'FontSize',14,'FontWeight','bold');
end

img=frame_curr;
%I=imcrop(img);
I=img;
I1=rgb2gray(I);
I2=adaptivethreshold(I1,15,0.02,0);
%I2=I1 <= 110;
D=bwdist(I2);
imshow(D)
L=watershed(-D);
imshow(L)
w=L==0;
figure, imshow(w);  
g2= I2 & ~w;
figure, imshow(g2);
