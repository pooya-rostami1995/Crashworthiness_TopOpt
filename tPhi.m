function [tmpPhi]=tPhi(xy,LSgridx,LSgridy,p)
st = xy(5);ct = sqrt(abs(1-st*st));
x1 = ct*(LSgridx - xy(1))+st*(LSgridy - xy(2));
y1 = -st*(LSgridx - xy(1))+ct*(LSgridy - xy(2));bb = xy(4);
tmpPhi = -((x1).^p/xy(3)^p+(y1).^p./bb.^p-1);
end