hairpins = "";
[rows,cols] = size(spacers);
for i=1:rows
    s.Sequence = convertStringsToChars(spacers(i,1));
    result = oligoprop(s);
    [m,n] = size(result.Hairpins);
    hairpins = hairpins + m;
end
fid = fopen('num_hairpins.txt','wt');
fprintf(fid, hairpins);
fclose(fid);

molweights = "";
for i=1:rows
    s.Sequence = convertStringsToChars(spacers(i,1));
    result = oligoprop(s);
    molweights = molweights + " " + result.MolWeight;
end
fid = fopen('mol_weights.txt','wt');
fprintf(fid, molweights);
fclose(fid);

deltaG = "";
for i=1:rows
    s.Sequence = convertStringsToChars(spacers(i,1));
    result = oligoprop(s);
    deltaG = deltaG + " " + result.Thermo(1,3);
end
fid = fopen('deltaG.txt','wt');
fprintf(fid, deltaG);
fclose(fid);


