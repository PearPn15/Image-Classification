function show_annotation(imgfile, annotation_file);
%% 
%% imgfile: string
%% annotation_file: string
%%
%% written by by Fei-Fei Li - November 2004
%%

IMTYPE = 'jpg'; 
GUIDELINE_MODE = 1; 
%% Parameters 
%label_abbrev = {'LE', 'RE', 'LN', 'NB', 'RN', 'LM', 'RM'}; 
LARGEFONT = 28; 
MEDFONT = 18; 
BIG_WINDOW = get(0,'ScreenSize'); 
SMALL_WINDOW = [100 100 512 480]; 
 
%% load the annotated data
load(annotation_file, 'box_coord', 'obj_contour');

%% Read and display image 
ima = imread(imgfile); 
ff=figure(1); clf; imagesc(ima); axis image; axis ij; hold on;
% black and white images
if length(size(ima))<3
   colormap(gray);
end
set(ff,'Position',SMALL_WINDOW); 
   
%% show box
box_handle = rectangle('position', [box_coord(3), box_coord(1), box_coord(4)-box_coord(3), box_coord(2)-box_coord(1)]);
set(box_handle, 'edgecolor','y', 'linewidth',5);

%% show contour
for cc = 1:size(obj_contour,2)
   if cc < size(obj_contour,2)
      plot([obj_contour(1,cc), obj_contour(1,cc+1)]+box_coord(3), [obj_contour(2,cc), obj_contour(2,cc+1)]+box_coord(1), 'r','linewidth',4);
   else
      plot([obj_contour(1,cc), obj_contour(1,1)]+box_coord(3), [obj_contour(2,cc), obj_contour(2,1)]+box_coord(1), 'r','linewidth',4);
   end
end

title(imgfile);
 
