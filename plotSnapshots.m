function plotSnapshots(filename, outfilename)

a = recolor(load(filename));
figure(1);
imshow(a,[]);
imwrite(a,outfilename);
