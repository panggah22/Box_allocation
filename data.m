function [box,container,c] = data(x)
switch x
    case 0
        box = [
            1   1   1 
            1   1   1   
            2   2   1]; 
        
        container = [
            2   2   2   5
            4   4   4   7];
    case 1
        box = [
            1 2 1
            1 2 1
            1 2 1
            2 2 2
            2 2 2
            2 2 2
            2 3 2
            2 3 2
            2 3 2
            2 4 2
            2 4 2
            2 4 2];
        
        container = [
            4 5 4 8
            4 5 4 8
            4 5 4 8
            4 6 4 10
            6 6 6 25];
    case 2
        box = [
            2 2 2
            2 2 2
            2 2 2
            2 2 2
            2 2 3
            2 2 3
            2 2 3
            2 2 3
            2 2 3
            2 2 3
            3 3 1
            3 3 1
            1 2 5];
        
        container = [
            3 3 7 80
            3 3 7 80
            4 4 7 110
            4 4 7 110];
            
end
box = box';
c = container(:,4)';
container = container(:,1:3)';
end