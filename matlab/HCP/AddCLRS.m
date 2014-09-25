function AddCLRS(infile1, infile2, outfile)

fid1 = fopen(infile1, 'r');
fid2 = fopen(infile2, 'r');
fidout = fopen(outfile, 'w');

clrs1 = cell2mat(textscan(fid1, '%f %f %f\n'));
clrs2 = cell2mat(textscan(fid2, '%f %f %f\n'));
fclose(fid1);
fclose(fid2);

clrs_out = (clrs1 + clrs2)/2;
clrs_out(sum(clrs1,2)==3,:) = clrs2(sum(clrs1,2)==3,:);
clrs_out(sum(clrs2,2)==3,:) = clrs1(sum(clrs2,2)==3,:);

for i = 1:size(clrs_out,1)
    fprintf(fidout, '%f %f %f\n', clrs_out(i,1), clrs_out(i,2), clrs_out(i,3));
end
fclose(fidout);
end