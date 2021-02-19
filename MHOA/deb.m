function [xb] = deb(xn, xo, D)

    if(xn(1,D+2) == xo(1,D+2))
        if(isnan(xn(1,D+1)) && isnan(xo(1,D+1)))
            if(rand() < 0.5)
                xb(1,:) = xn(1,:);
            else
                xb(1,:) = xo(1,:);
            end
        elseif(isnan(xn(1,D+1)))
            xb(1,:) = xo(1,:);
        elseif(isnan(xo(1,D+1)))
            xb(1,:) = xn(1,:);
        elseif(xn(1,D+1) < xo(1,D+1))
            xb(1,:) = xn(1,:);
        else
            xb(1,:) = xo(1,:);
        end
    elseif(xn(1,D+2) < xo(1,D+2))
        if(isnan(xn(1,D+1)))
            xb(1,:) = xo(1,:);
        else
            xb(1,:) = xn(1,:);
        end
    else
        if(isnan(xo(1,D+1)))
            xb(1,:) = xn(1,:);
        else
            xb(1,:) = xo(1,:);
        end
    end
end

% function xb = debo(xn, xo, D)
% 
%     if(xn(1,D+2)==0 && xo(1,D+2)~=0)
%         xb(1,:) = xn(1,:);
%     elseif(xn(D+2)~=0 && xo(1,D+2)==0)
%         xb(1,:) = xo(1,:);
%     elseif (xn(D+2)==0 && xo(1,D+2)==0)
%         if(xn(1,D+1) < xo(1,D+1))
%             xb(1,:) = xn(1,:);
%         else
%             xb(1,:) = xo(1,:);
%         end
%     else
%         if(xn(1,D+2) < xo(1,D+2))
%             xb(1,:) = xn(1,:);
%         else
%             xb(1,:) = xo(1,:);
%         end    
%     end
% end