

% This MATLAB SCRIPT reads the workspace generated in step1-script,and
% converts the raw data into actual platelet locations
% which is stored in COMdata.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
load('example_step1_workspace_RedPlt.mat');

List = cell(1, depth);
MajorList = [];
num = 0;
ind = 1;

FindNewObjects
   
for i = 1:length(w{1})

    List{1} = w{1}(i,2);

    for j = 1:length(w)

        ind = ind + 1;


        if j == 1

            List{ind} = w{j}(i,3);
            newName = w{j}(i,3);
%                List = [List, {newName}];

        else

%                index = find(w{j}(:,2) == newName);
%                [~, index, ~] = intersect(w{j}(:,2), newName) ;

            for loopIndex = 1:length(newName)

                if loopIndex == 1
                    index = [];
                end

                N = find(w{j}(:,2) == newName(loopIndex));
                index = [index; N];

            end

            if isempty(index)

%                    List = padarray(List, [0, 10-length(List)], 'post');
                emptyNode = cellfun('isempty', List);
                List(emptyNode) = {[0]};
                break;

            end

            newName = w{j}(index,3);
            List{ind} = w{j}(index,3);
%                List = [List, {newName}];

        end

    end

    MajorList = [MajorList; List];
    List = cell(1, depth);
    num = num + 1;
    ind = 1;
    
end

List = cell(1, depth);
SecondMajorList = [];
ind = 1;

for i = 1:length(Levels(1,:))
   
    for j = 1:sum(Levels(:,i)~=0)
        
%         List{1} = Levels(i,j);
        
        for k = (i + 1):length(w)
            
%             ind = ind + 1;
            
            if k == (i + 1)
                
                List(1:k-1) = {0};
                List{k} = Levels(j,i);    
                ind = k + 1;
                temp = w{k}(:,3);
                List{ind} = temp( (w{k}(:,2) == Levels(j,i)) );
                newName = List{ind};
            
            else
                
                ind = ind + 1;
            
                for loopIndex = 1:length(newName)

                    if loopIndex == 1
                        index = [];
                    end

                    N = find(w{k}(:,2) == newName(loopIndex));
                    index = [index; N];

                end
 
                if isempty(index)

                    emptyNode = cellfun('isempty', List);
                    List(emptyNode) = {[0]};
                    break;

                end

                newName = w{k}(index,3);
                List{ind} = w{k}(index,3);
            
            end
    
        end
        
    SecondMajorList = [SecondMajorList; List];
    List = cell(1, depth);
    num = num + 1;
    ind = 1;        
        
    end
    
end

MajorList = [MajorList; SecondMajorList];

%% this session only for plt, comment out this session for fibrin 

heightLimit =  22;  % use to seperate connected plt cluster
filamentCutoff = 3;

for i = 1:length(MajorList)
    count = 0;
    for j = 1:length(MajorList(1,:))
              
        if MajorList{i,j}  ~= 0
            count = count + 1;
        elseif count ~= 0
            break;
        else    
            continue;
        end
        
        if count > heightLimit
            count = 1;
            MajorList(end+1,:) = repmat({0}, 1, length(MajorList(1,:)));
            MajorList(end, (j-heightLimit):(j-1)) = MajorList(i, (j-heightLimit):(j-1));
            MajorList(i, (j-heightLimit):(j-1)) = repmat({0}, 1, heightLimit);
        end
        
    end
            
end

%% 
objData = ExtractFibrinNetworkData(MajorList, x, y, numObjs);
indicator = 1;
COMdata = zeros(length(objData(:,1,1)), 4);
filamentCutoff = 11;

for i = 1:length(objData(:,1,1))
    
    aveX = 0;
    aveY = 0;
    aveZ = 0;
    avePixels = 0;
    
    for j = 1:length(objData(1,:,1))
        
        if ~isempty(objData{i, j, 1})
            
            indicator = 1;
            saveJ = j;
            
            while indicator == 1
               
                for k = 1:length(objData{i,j,1})
                    
                    aveX = aveX + objData{i,j,1}(k).*objData{i,j,3}(k);
                    aveY = aveY + objData{i,j,2}(k).*objData{i,j,3}(k);
                    avePixels = avePixels + objData{i,j,3}(k);
                    aveZ = aveZ + j.*objData{i,j,3}(k);                   
                    
                end
                
                j = j + 1;
                
                if j > length(objData(1,:,1))
                    
                    indicator = -1;

                    height = j - saveJ;
                    
                    if height < filamentCutoff
                        
                        avePixels = 0;
                        
                    end
                    
                elseif isempty(objData{i,j,1})
                    
                    height = j - saveJ;
                    
                    if height < filamentCutoff
                        
                        avePixels = 0;
                        
                    end
                    
                    break;
                    
                end
                
            end
            
            break;
            
        end
               
    end
    
% if i == 4
%     
%     omega =1;
%     
% end

    COMdata(i,1) = aveX./avePixels;
    COMdata(i,2) = aveY./avePixels;
    COMdata(i,3) = 1.*aveZ./avePixels;
    COMdata(i,4) = avePixels;
    
    
end

radii = ((3)./(4.*pi)).*COMdata(:,4).^(1./3);
figure;
[x,y,z] = sphere;
FinalData = [];
COMcopy = COMdata;

% for i = 1:length(COMdata(:,1))
%  
%     newValues = COMdata(i,:);
%     
%     if (COMdata(i, 4) == 0)
%         continue;
%     end    
%         
%     for j = 1:length(COMdata(:,1))
%          
%         if (COMdata(j, 4) == 0) || j == i
%             continue;
%         end    
%         
%         distance = sqrt ( ( abs(COMdata(i,1) - COMdata(j,1) ) ).^2 ...
%                         + ( abs(COMdata(i,2) - COMdata(j,2) ) ).^2 ...
%                         + ( abs(COMdata(i,3) - COMdata(j,3) ) ).^2 );
%                     
%         if ( radii(i) + radii(j) ) <= distance
%             
%             continue;
%             
%         else 
% 
%             newPixels = newValues(1,4) + COMdata(j, 4);
%             aveX =  ( newValues(1,1).*newValues(1,4) + COMdata(j,1).*COMdata(j,4) )./ ( newPixels );
%             aveY =  ( newValues(1,2).*newValues(1,4) + COMdata(j,2).*COMdata(j,4) )./ ( newPixels );
%             aveZ =  ( newValues(1,3).*newValues(1,4) + COMdata(j,3).*COMdata(j,4) )./ ( newPixels );
%             newValues = [ aveX, aveY, aveZ, newPixels ];
%             COMdata(j,4) = 0;
%             radii(j) = 0;
%             
%         end
%                                                                                                                                                                                     
%     end
%     
%     COMdata(i,:) = newValues;
%     radii(i) = ((3.*pi)./4).*COMdata(i,4).^(1./3);
%     
% end
% 
% COMdata = COMdata(COMdata(:,4) ~= 0, :);
% radii = ((3.*pi)./4).*COMdata(:,4).^(1./3);

% for i = 1:length(COMdata(:,1))
% 
%     
%     
% end

for i = 1:length(COMdata(:,1))
    
    hold on
    surf((x.*radii(i) + COMdata(i,1)), (y.*radii(i) + COMdata(i,2)), (z.*radii(i) + COMdata(i,3)))
    
end

%         scatter3(x(w{i}(j,2)), y(w{i}(j,2)), i.*5, 'b', 'filled');
%         hold on
%         scatter3(x(w{i}(j,3)), y(w{i}(j,3)), (i + 1).*5, 'b', 'filled');
        
%         if j == 1
%             
%             old = w{j}(:,2);
%             new = w{j}(:,3);
%             
%         else
%            
%             nextLevel = w{j}(:,2);
%             common = intersect(old, nextLevel);
%             not_common = setdiff(old, nextLevel);
%             old = w{j}(:,2);
%             new = w{j}(:,3);
%             
%         end
            

    
%end