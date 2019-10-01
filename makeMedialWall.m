%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the index of medial wall for all parcellations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hem = {'L','R'};
maps={'Glasser_2016','Yeo_JNeurophysiol11_17Networks','Yeo_JNeurophysiol11_7Networks','Power2011','Yeo_CerCor2015_12Comp','Desikan','Dextrieux'};

for h=1:length(Hem)
    mwIdx = [];
    for i=1:length(maps)
        A=gifti(sprintf('%s.32k.%s.label.gii',maps{i}, Hem{h}));
        mwIdx = union(find(A.cdata(:,1)==0), mwIdx);
    end
    outfile = sprintf('medialWallIndex_%s.mat', Hem{h});
    save(outfile, 'mwIdx')
end