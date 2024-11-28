files = dir('*.bval');
for ij=1:length(files)
    bv = load(files(ij).name);
    bv(bv<=1) = 0;
    save(files(ij).name,'bv','-ascii');
end
