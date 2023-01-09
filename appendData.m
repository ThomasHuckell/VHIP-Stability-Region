function [oldData, lastEntry] = appendData(oldData,newData,startIndex)
% function that dynamically adds ode45 data, if startIndex is provided it will
% manually add to the data at the starting index

% get current size of data
dim = size(oldData);

% get new data size
newDim = size(newData);


switch nargin
    
    case 3
        
        if startIndex == 0
            oldData(startIndex + 1:startIndex + newDim(1),:) = newData(1:end,:);
        else
            
            oldData(startIndex + 1:startIndex + newDim(1)-1,:) = newData(2:end,:);
        
        end
        
    case 2
        
        % find last entry of the pre allocated data
        lastEntry = find(~ismember(oldData,zeros(1,dim(2)),'rows'),1,'last');
        
        % if first time adding data this is needed
        if isempty(lastEntry)
            
            oldData(1: newDim(1),:) = newData(1:end,:);
            
            lastEntry = 0;
        else
            % check if more space is need for data
            if lastEntry + newDim(1) > dim(1)
                
                % preallocate more space
                oldData(dim(1) + 1000, dim(2)) = 0;
                
            end
            
            % add new ode data to simdata
            oldData(lastEntry + 1:lastEntry + newDim(1)-1,:) = newData(2:end,:);
            
            
        end
        
        
        
end






end
