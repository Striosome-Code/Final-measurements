%{
Input: 
1) all output cropped images folders resulting from running the dots.m
script (dots to cropped images) on a folder of hemisphere images and
tetrode dot images.
2) excel file 'tetrodelist', filled in during dots.m  
3) new blank excel file 'measurements' -- must create this prior to running
4) relevant values, most at the beginning of the script, including: 
values of 0 or 1 to turn the distance measurement, selectivity measurement, and 
second desquig method on or off; 
diameters of squiggles and small regions to be removed; 
avg pixels per mm in the images; 
median gfp brightnesses for the gfp images; 
darkness threshold (search "dark") for finding ventricles

Ouput:
1) excel file 'measurements', filled in with relevant values for each
tetrode
2) 'output dot images' folder with MOR1 images showing a red dot at the
nearest striosome pixel to each tetrode
3) 'output BW images' folder showing desquigged striosome images with
ventricles highlighted in red 
%}

distance_switch = 1; %1 for on, 0 for off
selectivity_switch = 1; %" "
desquig2_switch = 0; %" "

pixels_per_mm = 1000; %fill in avg from images

small_region_diam = 10; %microns
squig_diam = 10; %" "

small_region_diam = small_region_diam*pixels_per_mm/1000; %pixels
squig_diam = squig_diam*pixels_per_mm/1000; %" "

data = xlsread('tetrodelist.xlsx');
ttnums = data(1,:);
ttregnums = data(2,:);

if selectivity_switch == 1
    gfpmedians = zeros(length(ttnums),7); %i.e. each row represents the medians
    %for a stack centered around a tetrode, going in order by tetrodes. 7 would
    %be the max num imgs in a stack. Values found using gfpmedians code beforehand
end

%---DENSITY (AND SELECTIVITY) MEASURE---
for c = 1:length(ttnums)
    files = dir(fullfile(['Output cropped images - tetrode ',num2str(ttnums(c))],'*.tif')); %tif files in folder 
    num_files = length(files); %number of files in folder 

    %cells to store different types of images and their names 
    mor1 = {};
    mor1names = {};
    strio = {};
    if selectivity_switch == 1
        gfp = {};
        gfpnames = {};
    end    
    ventricle_list = {};

    %go through images and store mor1, tt and their names
    for k = 1:num_files
        filename = files(k).name;
        [pathstr, name, ext] = fileparts(filename);    

        img = imread(fullfile(['Output cropped images - tetrode ',num2str(ttnums(c))],filename));
        img = rgb2gray(img); 
        
        if isempty(strfind(name, 'MOR1')) == 0
            mor1{end + 1} = img;    
            mor1names{1,end + 1} = name;
            
            dark = img < 10;    
            ventricles = bwareaopen(dark, round((small_region_diam/2)^2*pi));   
            ventricle_list{1,end + 1} = ventricles;

            img_not_ventricle = img(~ventricles);
            level = graythresh(img_not_ventricle);
            striosomes = im2bw(img, level);
            strio{end + 1} = striosomes;
            
        elseif isempty(strfind(name, 'GFP')) == 0
            if selectivity_switch == 1
                gfp{end + 1} = img;        
                gfpnames{1,end + 1} = name;  
            end
        end
    end;

    mor1nums = zeros(1,length(mor1)); %to store region numbers in MOR1 imgs
    if selectivity_switch == 1
      gfpnums = zeros(1,length(gfp));
    end

    %parse first 2 digits from image name string into region number for each
    %image type
    for k = 1:length(mor1nums)
        mor1nums(k) = str2num(mor1names{1,k}(1:2));
    end;
    
    if selectivity_switch == 1
        for k = 1:length(gfpnums)
            gfpnums(k) = str2num(gfpnames{1,k}(1:2));
        end
    end

    %find tetrode region number and corresponding MOR1 img
    ttregnum = ttregnums(c);
    mid_index = (mor1nums == ttregnum);
    mid_index = find(mid_index);
    mid_img = mor1{1,mid_index};

    [y,x] = size(mid_img); %flipped x and y because rows = y, cols = x
    [tetrode] = [round(x/2), round(y/2)]; %given it's at exact center of middle
    %image

    distances = zeros(1,floor(length(mor1)/2) + 1);
    for a = 1:floor(length(mor1)/2) + 1
        distances(a + 1) = 30*a; %in microns; fill in these values
    end;

    areas = zeros(1,length(distances) - 1); %to store total area across layers determined 
    %by each distance value
    num_strio_px = zeros(1,length(distances) - 1); %to store total number of striosome
    %pixels in layers determined by each distance value
    ratios = zeros(1,length(distances) - 1); %to store ratios of number of striosome
    %pixels to total area of all regions determined by each distance value

    for i = 1:length(distances) - 1
        %middle layer - consider as separate case
        r1 = pixels_per_mm/1000*distances(i); %converts distance to pixels
        r2 = pixels_per_mm/1000*distances(i + 1);
        middle_img = strio{mid_index};
        [x,y] = size(middle_img); 

        %if a pixel is within the circle, increment num_strio_px if its
        %brightness is above the threshold
        
        for k = 1:x
            for j = 1:y
                if (k - tetrode(1))^2 + (j - tetrode(2))^2 < r2^2 && (k - tetrode(1))^2 + (j - tetrode(2))^2 >= r1^2
                    num_strio_px(i) = num_strio_px(i) + (middle_img(k,j) == 1);
                end;
            end;
        end;

        areas(i) = areas(i) + pi*r2^2 - pi*r1^2;        


        if length(strio) > 1
            layers = floor(distances(i)/30); %how many layers above and below center 
            %to consider for area calculations
            for m = 1:layers                
                check1 = 0;
                check2 = 0;

                if any(mor1nums == ttregnum + m)
                    check1 = 1;
                    image1 = strio{mid_index + m};
                    [a,b] = size(image1);
                end;

                if any(mor1nums == ttregnum - m)                
                    check2 = 1;
                    image2 = strio{mid_index - m};
                    [cc,d] = size(image2); %(technically they should all be the same size
                    %but we don't want out-of-bounds errors)
                end;

                radius1 = (r1^2 - ((30*pixels_per_mm/1000)*m)^2)^(1/2); %pythagorean theorem
                radius2 = (r2^2 - ((30*pixels_per_mm/1000)*m)^2)^(1/2);
                areas(i) = areas(i) + 2*pi*radius2^2 - 2*pi*radius1^2; %because for each value of layers,
                %there are 2 circles
                
                %loop for image1
                if check1 == 1
                    for k = 1:a
                        for j = 1:b
                            if (k - tetrode(1))^2 + (j - tetrode(2))^2 < radius2^2 && (k - tetrode(1))^2 + (j - tetrode(2))^2 >= radius1^2
                                num_strio_px(i) = num_strio_px(i) + (image1(k,j) == 1);
                            end;
                        end;
                    end;
                end;

                %loop for image 2
                if check2 == 1
                    for k = 1:cc
                        for j = 1:d
                            if (k - tetrode(1))^2 + (j - tetrode(2))^2 < radius2^2 && (k - tetrode(1))^2 + (j - tetrode(2))^2 >= radius1^2
                                num_strio_px(i) = num_strio_px(i) + (image2(k,j) == 1);
                            end;
                        end;
                    end;
                end;
            end;                   
        end;
        ratios(i) = num_strio_px(i)/areas(i); %calculate/store ratio
    end;

    display(ratios);

    %cell to store distances & ratios in 2 columns 
    A = cell(1,3);
    B = cell(length(distances),3);
    
    A{1,1} = 'Tetrode';
    A{1,2} = 'Distance';
    A{1,3} = 'Density Ratio';
    
    B{1,1} = ttnums(c); %tetrode number
    
    for t = 1:length(distances) - 1 
        B{t,2} = distances(t + 1);
        B{t,3} = ratios(t);
    end;

    filename = 'measurements.xlsx'; %must have made this file beforehand
    xlswrite(filename,A); %write distances and ratios to excel file
    xlswrite(filename,B,1,['A',num2str((c - 1)*3 + 2)]);

    length_cst = 60; % length constant in microns
    length_cst_px = pixels_per_mm/1000*length_cst; %length constant in pixels
    [s,t] = size(strio{1});

    influence = zeros(1,s*t*length(strio)+100); %matrix to store exponential decay
    %influence values

    for i = 1:length(strio)
        [x,y] = size(strio{i});
        display(i);
        for k = 1:x
            for j = 1:y
                side1 = (abs(mid_index - i)*pixels_per_mm/1000*30); %distance in z direction
                side2 = ((abs(k - tetrode(1)))^2 + (abs(j - tetrode(2)))^2)^(1/2); %planar distance from tetrode tip
                distance = (side1^2 + side2^2)^(1/2); 
                px_influence = strio{i}(k,j) * (1 - exp(-distance/length_cst_px));
                influence(i*k+j) = px_influence;
            end;
        end;
    end;

    total_influence = sum(influence); %total tetrode influence with exponential
    %decay 
    display(total_influence);

    %write total tetrode influence to excel file
    C = {'With Exponential Decay'};
    xlswrite(filename,C,1,'D1');
    xlswrite(filename,total_influence,1,['D',num2str((c - 1)*3 + 2)]);

clear influence %bc this takes a ton of memory

%---DISTANCE MEASURE---
    if distance_switch == 1

        mkdir(['Output dot images - tetrode ',num2str(ttnums(c))]);

        smoothed = cell(1,length(strio)); %to store imgs with small regions removed 
        %and desquiggled with desquig
        smoothed2 = cell(1,length(strio)); %to store imgs with small regions removed 
        %and desquiggled with *desquig2*

        for h = 1:length(smoothed)
            smoothed{h} = bwareaopen(strio{h},round((small_region_diam/2)^2*pi));
            f = desquiggeneral(smoothed{h}, squig_diam, small_region_diam, 1);
            smoothed{h} = f;
            
            if desquig2_switch == 1
                smoothed2{h} = smoothed{h};            
                g = desquiggeneral(smoothed2{h}, squig_diam, small_region_diam, 2);
                smoothed2{h} = g;
            end          
        end;
    
        %e.g. if 7 layers, distances would be -90, -60, -30, 0, 30, 60, 90
        distances2 = zeros(1,length(distances)*2 - 3);
        distance_sort = sort(distances);
        for i = 1:length(distances2)
            distances2(i) = 30*(i-1) - distance_sort(end-1); %subtract 2nd largest
            %elt in distance matrix
        end;

        results = {}; %will store distance and layer of each striosome pixel
        %results2 = {}; %for second desquig method

        %FIRST DESQUIG METHOD
        for i = 1:length(distances2) %loops through images
            if any(mor1nums == ttregnum + distances2(i)/30)
                %img = smoothed{i};
                img = strio{i};
                [f,g] = size(img);
                for k = 1:f %loops through x coordinates
                    for j = 1:g %loops through y coordinates
                        if img(k,j) == 1 %check if above 
                            %striosome brightness threshold
                            side1 = pixels_per_mm/1000*abs(distances2(i)); %distance in z-dir
                            side2 = ((k - tetrode(1))^2 + (j - tetrode(2))^2)^(1/2); 
                            %distance in xy plane                
                            distance = (side1^2 + side2^2)^(1/2); %total distance
                            distance_mcr = distance/(pixels_per_mm/1000); %total distance 
                            locx = j;
                            locy = k;
                            %in microns
                            results{end + 1} = [i, distance_mcr, locx, locy]; %stores layer and distance
                        end;
                    end;
                end;
            end;
        end;

        %SECOND DESQUIG METHOD
        if desquig2_switch == 1
            for i = 1:length(distances2) %loops through images
                if any(mor1nums == mid_index + distances(i)/30)
                    img2 = smoothed2{i};
                    [f,g] = size(img2);
                    for k = 1:f %loops through x coordinates
                        for j = 1:g %loops through y coordinates
                            if img2(k,j) == 1 %check if above 
                                %striosome brightness threshold
                                side1 = pixels_per_mm/1000*abs(distances(i)); %distance in z-dir
                                side2 = ((k - tetrode(1))^2 + (j - tetrode(2))^2)^(1/2); 
                                %distance in xy plane                
                                distance = (side1^2 + side2^2)^(1/2); %total distance
                                distance_mcr = distance/(pixels_per_mm/1000); %total distance 
                                %in microns 
                                locx = j;
                                locy = k;
                                results2{end + 1} = [i, distance_mcr, locx, locy]; %stores layer and distance
                            end;
                        end;
                    end;
                end;
            end;
        end

        minimum = results{1,1}(1,2); %initializes minimum distance
        min_layer = results{1,1}(1,1); %initializes layer of minimum distance

        %DESQUIG2
        if desquig2_switch == 1
            min2 = results2{1,1}(1,2); %initializes minimum distance
            min_layer2 = results2{1,1}(1,1); %initializes layer of minimum distance
        end

        %loops through all striosome pixels and finds the minimum distance, along
        %with the layer in which it occurs and coordinates
        for h = 1:length(results)
            if results{1,h}(1,2) < minimum
                minimum = results{1,h}(1,2);
                min_layer = results{1,h}(1,1);
                min_locx = results{1,h}(1,3);
                min_locy = results{1,h}(1,4);
            end;
        end;

        %DESQUIG2
        if desquig2_switch == 1
            for h = 1:length(results2)
                if results2{1,h}(1,2) < min2
                    min2 = results2{1,h}(1,2);
                    min_layer2 = results2{1,h}(1,1);
                    min_locx2 = results2{1,h}(1,3);
                    min_locy2 = results2{1,h}(1,4);
                end;
            end;
        end

        display(minimum); %shows minimum distance from tetrode tip to striosome
        min_layer = mor1nums(min_layer);
        display(min_layer); %shows the number of the layer in which the minimum occurs

        %write min dist and min layer to excel file
        C = {'Minimum distance','Layer in stack'}; 
        D = {minimum, min_layer};
        xlswrite(filename,C,1,'E1');
        xlswrite(filename,D,1,['E',num2str((c - 1)*3 + 2)]);

        %find image containing minimum distance point
        min_ind = (mor1nums == min_layer);
        min_ind = find(min_ind);
        min_layer_img = mor1{min_ind};

        %plot minimum distance layer with red circle 10px diameter at point
        iptsetpref('ImshowBorder','tight');
        figure; imshow(min_layer_img);
        hold on 
        filledCircle([min_locx,min_locy],5,100,'r');
        h = gcf;
        saveas(h,['Output dot images - tetrode ',num2str(ttnums(c)),'/',mor1names{min_ind},'_dot.tif'],'tif'); %save 
        %image with dot in output folder
        hold off

        %DESQUIG2
        if desquig2_switch == 1
            display(min2); %shows minimum distance from tetrode tip to striosome
            min_layer2 = mor1nums(min_layer2);
            display(min_layer2); %shows the number of the layer in which the minimum occurs

            %write min dist and min layer to excel file
            C = {'Minimum distance - desquig2','Layer in stack - desquig2'; min2, min_layer2};
            xlswrite(filename,C,1,'F1');

            %find image containing minimum distance point
            min_ind2 = (mor1nums == min_layer2);
            min_ind2 = find(min_ind2);
            min_layer_img2 = mor1{min_ind2};

            %plot minimum distance layer with red circle 10px diameter at point
            iptsetpref('ImshowBorder','tight');
            figure; imshow(min_layer_img2);
            hold on 
            filledCircle([min_locx2,min_locy2],5,100,'r');
            h = gcf;
            saveas(h,['Output dot images - tetrode ',num2str(ttnums(c)),'/',mor1names{min_ind},'_dot_2.tif'],'tif'); %save 
            %image with dot in output folder
            hold off
        end
    end   

%---SELECTIVITY---
    %IMAGES
    if distance_switch ~= 1
        smoothed = cell(1,length(strio)); 
    end
    if selectivity_switch == 1
        mkdir(['Output BW images - tetrode ',num2str(ttnums(c))]);
        for k = 1:length(strio)
            if distance_switch == 1 %desquigging already done in distance section
                img = smoothed{k};
            else %desquig it now
                img = bwareaopen(strio{h},round((small_region_diam/2)^2*pi));
                f = desquiggeneral(img, squig_diam, small_region_diam, 1);
                smoothed{k} = f;                
            end
               [a,b] = size(img);
               rgbimg = repmat(double(img),[1 1 3]); %make RGB img
               current_ventricles = ventricle_list{1,k};  
               %loop through x,y coords
               for j = 1:b
                   for i = 1:a
                       if current_ventricles(j,i) == 1 %if current px is ventricle
                          rgbimg(j,i,:) = [1 0 0]; %color red
                       end
                   end
               end
               imwrite(rgbimg,['Output BW images - tetrode ',num2str(ttnums(c)),'/',mor1names{k},'_bw.tif'],'tif'); %save desquig BW img with red ventricles
        end     
        
        %MEASUREMENT
        striobrt = [];
        matrixbrt = [];

        current_gfpmedians = gfpmedians(c,:);
        
        middle_img = smoothed{mid_index};
        [x,y] = size(middle_img);
        mid_gfp = gfp{mid_index};
        midmedian = current_gfpmedians(mid_index);
        
        for i = 1:length(distances)
            %middle layer - consider as separate case
            r = pixels_per_mm/1000*distances(i); %converts distance to pixels
            middle_img = strio{mid_index} & ~ventricle_list{mid_index};
            [x,y] = size(middle_img); 

            %if a pixel is within the circle, increment num_strio_px if its
            %brightness is above the threshold
            for k = 1:x
                for j = 1:y
                    if (k - tetrode(1))^2 + (j - tetrode(2))^2 <= r^2
                        if middle_img(k,j) == 1
                            striobrt = [striobrt mid_gfp(k,j) - midmedian];
                        else
                            matrixbrt = [matrixbrt mid_gfp(k,j) - midmedian];
                        end
                    end;
                end;
            end;

            if length(strio) > 1
                layers = floor(distances(i)/30) - 1; %how many layers above and below center 
                %to consider for area calculations
                for m = 1:layers
                    check1 = 0;
                    check2 = 0;

                    if any(mor1nums == ttregnum + m) && any(gfpnums == ttregnum + m)
                        check1 = 1;
                        image1 = strio{mid_index + m} & ~ventricle_list{mid_index + m};
                        gfp1 = gfp{mid_index + m};
                        median1 = current_gfpmedians(mid_index + m);
                        [a,b] = size(image1);
                    end;
                    
                    if any(mor1nums == ttregnum - m) && any(gfpnums == ttregnum - m)              
                        check2 = 1;
                        image2 = strio{mid_index - m} & ~ventricle_list{mid_index - m};
                        gfp2 = gfp{mid_index - m};
                        median2 = current_gfpmedians(mid_index - m);
                        [cc,d] = size(image2); 
                    end;

                    radius = (r^2 - ((30*pixels_per_mm/1000)*m)^2)^(1/2); %pythagorean theorem
                   
                    %loop for image1
                    if check1 == 1
                        for k = 1:a
                            for j = 1:b
                                if (k - tetrode(1))^2 + (j - tetrode(2))^2 <= radius^2
                                    if image1(k,j) == 1
                                        striobrt = [striobrt (gfp1(k,j) - median1)];
                                    else
                                        matrixbrt = [matrixbrt (gfp1(k,j) - median1)];
                                    end
                                end;
                            end;
                        end;
                    end;

                    %loop for image 2
                    if check2 == 1
                        for k = 1:cc
                            for j = 1:d
                                if (k - tetrode(1))^2 + (j - tetrode(2))^2 <= radius^2
                                    if image2(k,j) == 1
                                        striobrt = [striobrt (gfp2(k,j) - median2)];
                                    else
                                        matrixbrt = [matrixbrt (gfp2(k,j) - median2)];
                                    end
                                end;
                            end;
                        end;
                    end;

                end;        
            end;
        end;
        
        %calculate means
        meanstrio = mean(striobrt(:));
        meanmatrix = mean(matrixbrt(:));

        %calculate ratio
        ratio = meanstrio/meanmatrix;
        
        %save ratio
        C = {'Selectivity Ratio'};
        xlswrite(filename,C,1,'G1');
        xlswrite(filename,ratio,1,['G',num2str((c - 1)*3 + 2)]);        
    end
end
