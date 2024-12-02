function cleanOutput(SITES)
    for s=1:length(SITES)
        for i=1:SITES{s}.ImageNum       
            system(['rm -r ',  SITES{s}.InputImages{i}.ImageDirectory, '/FSL/']);
        end
    end
end




