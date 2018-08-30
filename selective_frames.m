clc;
display('This project is about understanding the peristalsis motion of an earthworm');
%%To load the .data and .names files
display('Please browse through and select the right path that contains the video file');
path_folder=uigetdir('C:\', 'Select path to the data');
video_file_path= dir(fullfile(path_folder,'*.mp4'));
video_file_name=strsplit(video_file_path.name,'.');
video_name=video_file_name{1};
cd(path_folder);

ext=video_file_name{1,2};
movie_path_full=strcat(path_folder, '\', video_name, '.', ext);
video_object=VideoReader(movie_path_full);
v_h=video_object.Height;
v_w=video_object.Width;
frames_current=0;

frames_folder=strcat(path_folder,'\frames');
load ref_img;
load tail_ref;
%Testing only specific frames
framecount=0;
%num_frames=(5540-5500)+1;
%frame_curr1 = read(video_object,[5500 5540]);
num_frames=(6000-5500)+1;
frame_curr1 = read(video_object,[5500  6000]);
bound_cell={};
head_points_cell={};
tail_points_cell={};
clit_corn_cell={};
skelt_cell={};
bound_img_cell={};
L_cell={};
numm=0;

while numm<=(num_frames-1)
%try
framecount = framecount+1;
frame_curr = frame_curr1(:,:,:,numm+1);
numm=numm+25;
%tic
img=frame_curr;
I1=img;
C=I1;
I1=rgb2gray(I1);
I1=imtophat(I1,strel('disk',10));
I2=imadjust(I1);
%level=graythresh(I2);
%BW=im2bw(I2, (0.8*level));
BW=adaptivethreshold(I2,15,0.02,0);
BW2 = bwareaopen(BW, 20);
C(C<225)=0;
s=strel('disk',4,0);%Structuring element
D=~im2bw(C);%binary Image
% F=imerode(D,s);%Erode the image by structuring element
% filt=[1 1 1; 1 1 1; 1 1 1];
% erod_img=imerode(D,filt);
% bound_img=D-erod_img;
% bound_final=bwareaopen(bound_img,25);
% 
% diff_img=imadd(bound_final,BW2);
% bin_diff_img=im2bw(diff_img);
% Ibw=im2bw(I1);
% dif_t=(~Ibw)-bin_diff_img;

%Skeleton and boundary
area_f=imfill(D,'holes'); 
%boundary outline only
bound_out=bwmorph(area_f,'remove');

skel_f=bwmorph(area_f,'skel',Inf);

%Find the centroid points from final shape image :
s = regionprops(skel_f,'centroid');
centroids = cat(1, s.Centroid);

k=1;
for i=1:size(skel_f,1)
for j=1:size(skel_f,2)
if skel_f(i,j)==1
xdata_new(k)=i;
ydata_new(k)=j; k=k+1;
end
end
end

%smooth data
xo=xdata_new;
yo=ydata_new;
x=xdata_new;
y=ydata_new;
x=[x(end) x x(1)];
y=[y(end) y y(1)];
xs=zeros(1,length(xo));
ys=zeros(1,length(yo));
for i=2:length(x)-1
    ab=[x(i+1)-x(i-1); y(i+1)-y(i-1)];
    abu=ab/sqrt(ab(1)^2+ab(2)^2); 
    ac=[x(i)-x(i-1); y(i)-y(i-1)];
    d=dot(ac,abu)*abu+[x(i-1); y(i-1)];  
    [ex,ey]=midpoint(d(1),d(2),x(i),y(i));
    xs(i-1)=ex;
    ys(i-1)=ey;
end
%x and y coordinates of area_f
k=1;
for i=1:size(area_f,1)
for j=1:size(area_f,2)
if area_f(i,j)==1
xdata_area(k)=i;
ydata_area(k)=j; k=k+1;
end
end
end

%Getting the corner points from area_f (Looks good for this particular
%frame):
C = corner(area_f,'MinimumEigenvalue');

%Reducing the centroid points: 
s = regionprops(area_f,'centroid');
centroids = cat(1, s.Centroid);

k=1; 
for i=1:size(centroids,1)
%disp(area_f(round(centroids(i,2)),round(centroids(i,1))));
if(area_f(round(centroids(i,2)),round(centroids(i,1)))==1)
 centroids_final(k,1)=centroids(i,1);
centroids_final(k,2)=centroids(i,2); k=k+1;
end
end

%Find the head point using integral Image:
inttrial=integralImage(bound_out);
inttrial_bw=im2bw(inttrial); 
[rows,cols] = size(inttrial_bw);
for col = 1:cols
if size(unique(inttrial_bw(:,col)))==1
continue;
else
for row = 1:rows
if inttrial_bw(row,col)==1
    head_pt_x=col;
    head_pt_g=row;
    break;
end
end
end
break;
end
head_point_final=horzcat(head_pt_x,head_pt_g);

%Find the tail point:
c=normxcorr2(tail_ref, bound_out);
[ypeak, xpeak] = find(c==max(c(:)),1);
yoffSet = ypeak-size(tail_ref,1);
xoffSet = xpeak-size(tail_ref,2);
% hFig = figure;
% hAx  = axes;
% imshow(bound_out,'Parent', hAx);
% imrect(hAx, [xoffSet+1, yoffSet+1, size(tail_ref,2), size(tail_ref,1)]);
% title('Identified tail area by NCC');

k=1; [rows1,cols1] = size(tail_ref);
tail_pos_x=zeros;
tail_pos_y=zeros;

for row = 1:rows1
for col = 1:cols1
tail_pos_x(k)=(xoffSet+1)+row;
tail_pos_y(k)=(yoffSet+1)+col;
k=k+1;
end
end
tail_pos_x_unq=unique(tail_pos_x);
tail_pos_y_unq=unique(tail_pos_y); 
tail_final=ismember(C,horzcat(tail_pos_x_unq,tail_pos_y_unq));
k=1;
for i=1:size(C,1)
if (tail_final(i,1)==1 && tail_final(i,2)==1)
ind_x(k)=i;k=k+1;
end
end
tail_points_final=zeros;
for i=1:size(ind_x,2)
tail_points_final(i,1)=C(ind_x(i),1);
tail_points_final(i,2)=C(ind_x(i),2);
end

%get the clitellum co-ordinates by using SURF detection :
clit_bbox=SURF_BasedDetection(ref_img,frame_curr);
clit_pt_x=(min(clit_bbox(:,1)))+(max(clit_bbox(:,1))-min(clit_bbox(:,1)))/2;
clit_pt_y=(min(clit_bbox(:,2)))+(max(clit_bbox(:,2))-min(clit_bbox(:,2)))/2;
%figure, imshow(bound_out); hold on; plot(clit_pt_x,clit_pt_y,'r*');hold off

%Now plot all these points and the boundary :
[boundy,L] = bwboundaries(bound_out,'noholes');
%To save all important parameters for later calculation in cell:
bound_cell{framecount}=boundy;
head_points_cell{framecount}=head_point_final;
tail_points_cell{framecount}=tail_points_final;
clit_corn_cell{framecount}=horzcat(clit_pt_x,clit_pt_y);
skelt_cell{framecount}=skel_f;
bound_img_cell{framecount}=bound_out;
L_cell{framecount}=L;
imshow(imsharpen(imfuse(bound_out,skel_f)));
drawnow;
hold on;
for k = 1:length(boundy)
   boundary = boundy{k};
   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2);
   drawnow;
end
hold on; plot(head_point_final(:,1),head_point_final(:,2),'c*');
drawnow;
 hold on; plot(tail_points_final(:,1),tail_points_final(:,2),'c*');
 drawnow;
 hold on; plot(clit_pt_x,clit_pt_y,'r*');
 drawnow;
 hold on; 
 plot(centroids_final(:,1),centroids_final(:,2),'b*');
 drawnow;
 hold on;
%toc
% ME
 %   fprintf('Mostly the error would have been due to absence of clitellum: %s\n', ME.message);
 %   continue;  % Jump to next iteration of: for i
%end
end

%For plotting the area:
ll=zeros;
xxax=zeros;
for i=1:(framecount-1)
area_curr=regionprops(L_cell{1,i},'area'); area_cell=struct2cell(area_curr);
ll(i,1)=area_cell{1,1}; xxax(i,1)=i;
end;
plot(xxax,ll)

%For getting the length of the worm:
worm_length=zeros;
for i=1:(framecount-1)
skel_1=skelt_cell{1,i};
[row_1,col_1]=find(skel_1==1);
dist=0;
for j=1:(size(col_1,1)-1)
dist=dist+sqrt((row_1(j,1)-row_1(j+1,1))^2+(col_1(j,1)-col_1(j+1,1))^2);
end
worm_length(i)=dist;
end

%Head tail position variations :
hedd=zeros;
tl=zeros;
for i=1:(framecount-1)
hedd(i,1)=head_points_cell{1,i}(:,1); hedd(i,2)=head_points_cell{1,i}(:,2);
tl(i,1)=tail_points_cell{1,i}(1,1); tl(i,2)=tail_points_cell{1,i}(1,2);
end;
plot(hedd(:,1),hedd(:,2),'c*-');
hold on;
plot(tl(:,1),tl(:,2),'g*-');
