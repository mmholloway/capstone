function [peaks,loc,c, peakLocMatrix] = Find_Peaks(cells, loopnumber, size_limit, minValue)
   long = length(cells);
   high = length(cells(:,1));
   peakLocMatrix = zeros(long,high);
   cellAvg = zeros(long, high);

   cellAvghold = zeros(long, high);
   loc = zeros(length(cells)*3,2);
    for x=5:long-4
        for y=5:high-4
            if cells(x,y) >=minValue && cells(x,y)~=0
                cellAvg(x,y) = 1;
            end
        end
    end
    for loop=1:loopnumber
        for x=5:long-4
            for y=5:high-4
                cellAvghold(x,y) = cellAvg(x+1,y+1) + cellAvg(x+1,y) + cellAvg(x+1,y-1) + ...
                                   cellAvg(x,y+1) + cellAvg(x,y) + cellAvg(x,y-1) + ...
                                   cellAvg(x-1,y+1) + cellAvg(x-1,y) + cellAvg(x-1,y-1);
            end
        end
        cellAvg(:,:) = cellAvghold(:,:);
    end
    

   % rewrite for speed
   c=0;
   for x=10:long-9
       for y=10:high-9
            if cellAvg(x,y)>=cellAvg(x+1,y+1) && cellAvg(x,y)>=cellAvg(x+1,y) && cellAvg(x,y)>=cellAvg(x+1,y-1) && ...
               cellAvg(x,y)>=cellAvg(x,y+1) && cellAvg(x,y)>=cellAvg(x,y-1) && ...
               cellAvg(x,y)>=cellAvg(x-1,y+1) && cellAvg(x,y)>=cellAvg(x-1,y) && cellAvg(x,y)>=cellAvg(x-1,y-1) && ...
               cellAvg(x,y)>(size_limit*8^(loopnumber-1))
%                cellAvg(x,y)>=loopnumber
                c=c+1;
                loc(c,1)=x;
                loc(c,2)=y;
            end
       end
   end
   cmax = c;
   for original = 1:length(loc)
       if loc(original,1)>0
%            for check = 1:length(loc)
%                xFlag = false;
%                yFlag = false;
%                if check~=original
%                    for range = -round(size_limit/2):round(size_limit/2)
%                        if loc(original,1) == loc(check,1) + range
%                            xFlag = true;
%                        end
%                        if loc(original,2) == loc(check,2) + range
%                            yFlag = true;
%                        end
%                    end
%                    if xFlag && yFlag
%                        c=c-1;
%                        loc(check,1) = 0;
%                        loc(check,2) = 0;
%                    end
                   viewRange = cells(loc(original,1)-round(size_limit/2):loc(original,1)+round(size_limit/2),loc(original,2)-round(size_limit/2):loc(original,2)+round(size_limit/2));
                   for x = round(size_limit/2):length(viewRange)-round(size_limit/2)
                       for y = round(size_limit/2):length(viewRange(1,:))-round(size_limit/2)
                           if viewRange(x,y) > 0 && x ~= 1+size_limit/2
                               cells(loc(original,1)-round(size_limit/2)+x, loc(original,2)-round(size_limit/2)+y) = 0;
                           end
                       end
                   end
%                end
%            end
       end
   end
   for counter = 1:cmax
       for order = 2:length(loc)
           if loc(order,1)>loc(order-1,1)
               xholder = loc(order,1);
               loc(order,1) = loc(order-1,1);
               loc(order-1,1) = xholder;
               yholder = loc(order,2);
               loc(order,2) = loc(order-1,2);
               loc(order-1,2) = yholder;
           end
       end
   end
%    for x=1:length(cellAvg)
%        for y=1:length(cellAvg(1,:))
%            for i=1:c
%                if x==loc(c,1) && y==loc(c,2)
%                    peakLocMatrix(x,y) = 1;
%                end
%            end
%        end
%    end
   for i = 1:c
       peakLocMatrix(loc(i,1),loc(i,2)) = 1;
   end
   peaks = cellAvg;
end




