leftI = (imread("images\middlebury_left.png"));
rightI = (imread("images\middlebury_right.png"));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: to get reasonable runtimes I added code to downsample the image so
% that the algorithm can be quickly validated, otherwise the output would
% take a very long time
% 
% To replicate the results presented in the report for the fullsize image
% do the following:
%
% 1-Comment out the code the downsamples the images
% 2-In getDispMapLtoR and getDispMapRtoL change the disparity=0:127 to disparity=0:255
% 3-in fastBitwiseVoting change iter=0:6 to iter=0:7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Resize image for fast runtime
rightI=imresize(rightI, [375, 450], 'bilinear');
leftI=imresize(leftI, [375, 450], 'bilinear');

%number of iterations for the refinement process
iterations=2;


figure;
dispMapL=getDispMapLtoR(leftI,rightI,8);
imagesc(dispMapL)
title('dispmap L');
colormap("gray");

figure;
dispMapR=getDispMapRtoL(leftI,rightI,8);
imagesc(dispMapR);
title('dispmap R');
colormap("gray");

%This loop just computes and displays the results
for i=1:iterations
    figure
    [l,r] =disparityCrossCheck(dispMapL,dispMapR);
    dispMapL=l;
    dispMapR =r;
    imagesc(dispMapL)
    colormap("gray")
    titleText = sprintf('dispmap L after cross checking Iteration=%d', i);
    title(titleText);
    
    figure;
    imagesc(dispMapR)
    colormap("gray")
    titleText = sprintf('dispmap R after cross checking Iteration=%d', i);
    title(titleText);
    
    figure;
    dispMapL=fastBitwiseVoting(dispMapL,leftI);
    imagesc(dispMapL);
    titleText = sprintf('dispmap L after cross checking and bitwise voting Iteration=%d', i);
    title(titleText);
    colormap("gray")
    
    figure
    dispMapR = fastBitwiseVoting(dispMapR,rightI);
    imagesc(dispMapR);
    titleText = sprintf('dispmap R after cross checking and bitwise voting Iteration=%d', i);
    title(titleText);
    colormap("gray")
end


%This function calculates the disparity map by matching pixels from the
%left image to ones in the right image with a window of size w*2+1

%For faster run time with smaller images lower the range of disparities in
%the for loop
function dispMap = getDispMapLtoR(leftI,rightI,w)

    dispMap = zeros(size(leftI,1),size(leftI,2));    
    for r=1+w:1:size(leftI,1)-w
        r
        for c=size(leftI,2)-w:-1:1+w
            lastssd = inf;

            blockLeftR = double(leftI(r - w:r + w, c - w:c + w,1));
            blockLeftG = double(leftI(r - w:r + w, c - w:c + w,2));
            blockLeftB = double(leftI(r - w:r + w, c - w:c + w,3));

            for disparity=0:255
                if c-disparity<1+w
                    break;
                end
                
                blockRightR = double(rightI(r - w:r + w, c - w-disparity:c + w-disparity,1));
                blockRightG = double(rightI(r - w:r + w, c - w-disparity:c + w-disparity,2));
                blockRightB = double(rightI(r - w:r + w, c - w-disparity:c + w-disparity,3));
    
                ssd= sum(sum((blockLeftR-blockRightR).^2))+sum(sum((blockLeftG-blockRightG).^2))+sum(sum((blockLeftB-blockRightB).^2)) + (dispMap(r,c+1)-disparity)^2 + (dispMap(r-1,c)-disparity)^2;
                
                if ssd <lastssd
                    dispMap(r,c) = disparity;
                    lastssd=ssd;
                end         
            end
        end
    end
end

%Same as above but calculates the disparity map by matching pixels in the
%right image to ones in the left one
function dispMap = getDispMapRtoL(leftI,rightI,w)
  
    dispMap = zeros(size(leftI,1),size(leftI,2));        
    for r=1+w:1:size(leftI,1)-w
        r
        for c=size(leftI,2)-w:-1:1+w
            lastssd = inf;
            blockRightR = double(rightI(r - w:r + w, c - w:c + w,1));
            blockRightG = double(rightI(r - w:r + w, c - w:c + w,2));
            blockRightB = double(rightI(r - w:r + w, c - w:c + w,3));
            for disparity=0:255
                if c+disparity>size(leftI,2)-w
                    break;
                end
          
                blockLeftR = double(leftI(r - w:r + w, c - w+disparity:c + w+disparity,1));
                blockLeftG = double(leftI(r - w:r + w, c - w+disparity:c + w+disparity,2));
                blockLeftB = double(leftI(r - w:r + w, c - w+disparity:c + w+disparity,3));

                ssd= sum(sum((blockLeftR-blockRightR).^2))+sum(sum((blockLeftG-blockRightG).^2))+sum(sum((blockLeftB-blockRightB).^2)) + (dispMap(r,c+1)-disparity)^2 + (dispMap(r-1,c)-disparity)^2;

                if ssd <lastssd
                    dispMap(r,c) = disparity;
                    lastssd=ssd; 
                end 
            end
        end
    end
end

% Cross checks the disparities between the two disparity maps
% returns the two disparity maps with half occluded points assigned the
% value 0.
function [dispMap1,dispMap2] = disparityCrossCheck(dispMap1,dispMap2)
   
    for r=1:size(dispMap1,1)
        for c=1:size(dispMap1,2)
            if c-dispMap1(r,c)<1 || abs(dispMap1(r,c)- dispMap2(r,c-dispMap1(r,c))) ~=0
                dispMap1(r,c)=0;
                
            end
            if c+dispMap2(r,c)>size(dispMap1,2) || abs(dispMap2(r,c)- dispMap1(r,c+dispMap2(r,c))) ~=0
                dispMap2(r,c)=0;
                
            end
        end
    end

end

% Implements the process of fast bitwise voting given a disparity map and
% an image
function newDispMap = fastBitwiseVoting(dispMap,leftImage)
    newDispMap = dispMap;
    axesMatrix =getWindowAxes(leftImage);
    
    for row=1:size(dispMap,1)
        row
    for col =1:size(dispMap,2)
    axes = axesMatrix{row,col};
    horizontalSupportWindow = getHSW(row,col,leftImage,axes(4),axes(3));
    verticalSupportWindow = getVSW(row,col,leftImage,axes(2),axes(1));
    for iter=0:7
        counter =0;
        counterV=0;
        validEntryCounter=0;
        validEntryCounterV=0;
        mask = (2^iter);
        
        for i=1:size(horizontalSupportWindow,1)
            r=horizontalSupportWindow(i,1);
            for c=horizontalSupportWindow(i,2):1:horizontalSupportWindow(i,3)
                if int16(bitand(int16(dispMap(r,c)), int16(mask)))>0
                    counter =counter+1;
                end
                if dispMap(r,c)~=0
                    validEntryCounter = validEntryCounter+1;
                end
            end
        end
        
        for i=1:size(verticalSupportWindow,1)
            c=verticalSupportWindow(i,1);
            for r=verticalSupportWindow(i,3):1:verticalSupportWindow(i,2)
                if int16(bitand(int16(dispMap(r,c)), int16(mask)))>0
                    counterV=counterV+1;
                end
                
                if dispMap(r,c)~=0

                    validEntryCounterV = validEntryCounterV+1;
                end
            end
        end
        if 0.5*counter+0.5*counterV>0.5*((0.5*validEntryCounter)+(0.5*validEntryCounterV))
            newDispMap(row,col) =int16(bitor(int16(newDispMap(row,col)), int16(mask)));
        else   
            newDispMap(row,col) =int16(bitset(newDispMap(row,col),iter+1,0)); 
        end

    end
    end
    end

end



%creating Edge-Sensitive Local Support axes
% Threshold is defined by 'tau'
function supportWindowAxesMatrix = getWindowAxes(image)
    supportWindowAxesMatrix = cell(size(image,1),size(image,2));
    rows = size(image,1);
    cols = size(image,2);
    
    tau =30;
    for r=1:rows
        for c=1:cols

            quadruplet = [0,0,0,0];
            %finding the upward vectical bound
            for r_=r:-1:1
                rDiff  =abs(int16(image(r_,c,1))-int16(image(r,c,1)));
                gDiff  =abs(int16(image(r_,c,2))-int16(image(r,c,2)));
                bDiff  =abs(int16(image(r_,c,3))-int16(image(r,c,3)));
                maxDiff = max(max(rDiff,gDiff),bDiff);
                if maxDiff<= tau %&& r-r_<9
                    quadruplet(1,4)=r-r_;
                else
                    break;
                end
            end

            %finding the downward vertical bound
            for r_=r:1:rows
                rDiff  =abs(int16(image(r_,c,1))-int16(image(r,c,1)));
                gDiff  =abs(int16(image(r_,c,2))-int16(image(r,c,2)));
                bDiff  =abs(int16(image(r_,c,3))-int16(image(r,c,3)));
                maxDiff = max(max(rDiff,gDiff),bDiff);
                if maxDiff<= tau %&& r_-r<9
                    quadruplet(1,3)=r_-r;
                else
                    break;
                end
            end

             %finding the leftward horizontal bound
            for c_=c:-1:1
                rDiff  =abs(int16(image(r,c_,1))-int16(image(r,c,1)));
                gDiff  =abs(int16(image(r,c_,2))-int16(image(r,c,2)));
                bDiff  =abs(int16(image(r,c_,3))-int16(image(r,c,3)));
                maxDiff = max(max(rDiff,gDiff),bDiff);
                if maxDiff<= tau %&& c-c_<9
                    quadruplet(1,1)=c-c_;
                else
                    break;
                end
            end

            %finding the rightward horizontal bound
            for c_=c:1:cols
                rDiff  =abs(int16(image(r,c_,1))-int16(image(r,c,1)));
                gDiff  =abs(int16(image(r,c_,2))-int16(image(r,c,2)));
                bDiff  =abs(int16(image(r,c_,3))-int16(image(r,c,3)));
                maxDiff = max(max(rDiff,gDiff),bDiff);
                if maxDiff<= tau %&&c_-c<9
                    quadruplet(1,2)=c_-c;
                else
                    break;
                end
            end
            
            supportWindowAxesMatrix{r,c} = quadruplet;
        end
    end
end





%Gets the horizontal support window for a pixel given the image and the
%bounds for its vertical support axis
function horizontalSupportWindow = getHSW(r,c,image,Vplus,Vminus)
    horizontalSupportWindow = zeros(Vplus+Vminus+1,3);
    tau =30;
    startR = r-Vplus;

    r=startR;
    %Finding the leftward horizontal bound
    for v = 1:Vplus+Vminus+1
        horizontalSupportWindow(v,1) = r;
        for c_ = c:-1:1
            rDiff  =abs(int16(image(r,c_,1))-int16(image(r,c,1)));
            gDiff  =abs(int16(image(r,c_,2))-int16(image(r,c,2)));
            bDiff  =abs(int16(image(r,c_,3))-int16(image(r,c,3)));
            maxDiff = max(max(rDiff,gDiff),bDiff);
            if maxDiff<= tau %&& c-c_<9
                horizontalSupportWindow(v,2) = c_;
            else
                break;
            end
        end
        r=r+1;
    end
    
    r=startR;
    %Finding the rightward horizontal bound
    for v = 1:Vplus+Vminus+1
        horizontalSupportWindow(v,1) = r;
        for c_ = c:1:size(image,2)
            rDiff  =abs(int16(image(r,c_,1))-int16(image(r,c,1)));
            gDiff  =abs(int16(image(r,c_,2))-int16(image(r,c,2)));
            bDiff  =abs(int16(image(r,c_,3))-int16(image(r,c,3)));
            maxDiff = max(max(rDiff,gDiff),bDiff);
            if maxDiff<= tau %&& c_-c<9
                horizontalSupportWindow(v,3) = c_;
            else
                break;
            end
        end
        r=r+1;
    end
end

%Gets the vertical support window for a pixel given the image and the
%bounds for its horizontal support axis
function vecticalSupportWindow = getVSW(r,c,image,Hplus,Hminus)
    vecticalSupportWindow = zeros(Hplus+Hminus+1,3);
    tau =30;
    
    startC =c-Hminus;
    
    
    c=startC;
    %find upward vertical bound 
    for h = 1:Hplus+Hminus+1
        vecticalSupportWindow(h,1)=c;
        for r_ = r:-1:1
            rDiff  =abs(int16(image(r_,c,1))-int16(image(r,c,1)));
            gDiff  =abs(int16(image(r_,c,2))-int16(image(r,c,2)));
            bDiff  =abs(int16(image(r_,c,3))-int16(image(r,c,3)));
            maxDiff = max(max(rDiff,gDiff),bDiff);
            if maxDiff<= tau %&& r-r_<9
                vecticalSupportWindow(h,3) = r_;
            else
                break;
            end
        end
        c=c+1;
    end

    c=startC;
    %Find the downward vertical bound 
    for h = 1:Hplus+Hminus+1
        vecticalSupportWindow(h,1)=c;
        for r_ = r:1:size(image,1)
            rDiff  =abs(int16(image(r_,c,1))-int16(image(r,c,1)));
            gDiff  =abs(int16(image(r_,c,2))-int16(image(r,c,2)));
            bDiff  =abs(int16(image(r_,c,3))-int16(image(r,c,3)));
            maxDiff = max(max(rDiff,gDiff),bDiff);
            if maxDiff<= tau %&&r_-r<9
                vecticalSupportWindow(h,2) = r_;
            else
                break;
            end
        end
        c=c+1;
    end
end
