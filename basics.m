L = iread('Aloe\view1.png');
R = iread('Aloe\view5.png');
stdisp(L, R);

GT = iread('Aloe\disp1.png');
GT = GT / 2;
dmin = double(min(min(GT)));
dmax = double(max(max(GT)));
H = 3;
d = istereo(L, R, [dmin, dmax], H);

fileID = fopen('Aloe\dmin.txt');
offset = fscanf(fileID, '%d') / 2
d = d + offset;
GT = GT + offset;

figure
idisp(d, 'bar');
figure
stdisp(d, GT);