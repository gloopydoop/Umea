clear all
clc

rho = ones(5,9);
rho(:,1:3) = 0;
x = 5;
y = 3;
i = 9;
j = 2;

% find the smallest ball with (x,y) in it
r_first = ceil(norm([(x-i), (y-j)],2));
r = r_first - 1;
finished = false;

while finished == false
    % search area
    r = r+1;
    N = min(size(rho,1),j+r);
    S = max(1,j-r);
    E = min(size(rho,2),i+r);
    W = max(1,i-r);
    
    for m = W:E
        for n = S:N
            if norm([(n-j), (m-i)],2)<=r && rho(n,m) == 0
                finished = true;
            end
        end
    end
    
    % also break if we're at the wall
    if (N == size(rho,1) || S == 1 || E == size(rho,2) || W == 1)
        if finished == false
            r = r + 1;
        end
        finished = true;
    end
end

if r == r_first
    disp('failed_first')
else
    disp(num2str(r-1)) % it failed on the r'th. so the r
end
                
            
        
