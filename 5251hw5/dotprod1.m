function out=dotprod1(x,y);

col=size(x,2);
xs=x.*y;
out=cumsum(xs')';
out=out(:,col);
