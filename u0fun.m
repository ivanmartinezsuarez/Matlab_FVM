function [y] = u0fun(x)
%y=((-3<x)&&(x<3))*(sin(x))+0;
%y=log(abs(cos(1/x)));
%y=sin(cos(tan(abs(log(1./x.^3)))));
y=sin(x);
end

