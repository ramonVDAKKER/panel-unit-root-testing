function [yast,xast]=PUR_max(y)
%quadratic approximation of hight of top of 'function' y

[ymax,ypos]=max(y); 
xdelta=1;
if or(ypos==1,ypos==length(y))
%   fprintf('WARNING: max on boundary %10.0f\n',ypos)
   yast=ymax;
   xast=ypos; %on yscale 1…length(y)
else 
   F=y(ypos-1:ypos+1); 
   yast=F(2)-(F(1)-F(3))^2/8/(F(1)+F(3)-2*F(2)); 
   xast=ypos+xdelta/2*(F(1)-F(3))/(F(1)+F(3)-2*F(2)); %on yscale 1…length(y)
end;
