function objData = ExtractFibrinNetworkData(MajorList, x, y, numObjs)

objData = cell(length([MajorList{:,1}]), length([MajorList{1,:}]), 3);

for i = 1:length([MajorList{:,1}])
    
    for j = 1:length([MajorList{1,:}])
               
        if MajorList{i,j} ~= 0
                
            objData{i,j,1} = x(MajorList{i,j}, j);
            objData{i,j,2} = y(MajorList{i,j}, j);
            objData{i,j,3} = numObjs(MajorList{i,j}, j);
            
        else
            
            objData{i,j,1} = [];
            objData{i,j,2} = [];
            objData{i,j,3} = [];
            
%            break;
        
        end
            
    end
end