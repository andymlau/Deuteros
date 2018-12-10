function nummat = cell2num(cellofstringsofnum)
% Takes a cell of strings of num and unpacks to matrix of nums

nummat = [];
for i=1:length(cellofstringsofnum)
    line = cellofstringsofnum(i,:);
    
    new_line = [];
    for j=1:length(line)
        ele = str2num(cell2mat(line(j)));
        if isempty(ele)
            ele = 0;
        end
        new_line = [new_line; ele];
    end
    
    nummat = [nummat; new_line'];
end



