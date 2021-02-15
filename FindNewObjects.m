Levels = [];

for i = 1:(length(w(1,:)) - 1)
    
    if i == 119
        omega = 1;
    end
    
    currentLevel = w{i}(:,3);
    nextLevel = w{i+1}(:,2);
    not_common = [setdiff(nextLevel, currentLevel)];
    Levels(1:length(not_common), i) = not_common;
    
end