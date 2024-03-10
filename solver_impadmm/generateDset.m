% suppose volexs of size:   sizex * sizey  (suppose sizex = sizey)
function [Dset,K1,K2,K3] = generateDset(sizex,sizey)
if sizex == sizey
    Dset = cell(2*sizex-1,2);
    % D0,D-1,D1,D-2,D2,D-3,D3,...
    for ii = 1:sizex
        if ii == 1
            Dset{1,1} = 0;
            Dset{1,2} = [(1:sizex)',(1:sizex)'];
        else
            Dset{2*ii-2,1} = -(ii-1);
            Dset{2*ii-2,2} = [(ii:sizex)',(ii-1:sizex-1)'];
            Dset{2*ii - 1,1} = ii-1;
        %    Dset{2*ii - 1,2} = [(ii-1:sizex-1)',(ii:sizex)'];
            Dset{2*ii - 1,2} = [(1:sizex-ii+1)',(ii:sizex)'];
        end
    end
elseif sizex < sizey
    Dset = cell(sizex + sizey -1,2);
    % D0,D-1,D1,D-2,D2,D-3,D3,...Dsizex,Dsizex+1,...
    for ii = 1:sizey
        if ii == 1
            Dset{1,1} = 0;
            Dset{1,2} = [(1:sizex)',(1:sizex)'];
        elseif ii <= sizex
            Dset{2*ii-2,1} = -(ii-1);
            Dset{2*ii-2,2} = [(ii:sizex)',(1:sizex-ii+1)'];
            Dset{2*ii - 1,1} = ii-1;
            Dset{2*ii - 1,2} = [(1:min(sizey,sizex + ii-1)-ii+1)',(ii:min(sizey,sizex + ii-1))'];
        else
            Dset{ii + sizex - 1,1} = ii - 1 ;
            Dset{ii + sizex - 1,2} = [(1:min(sizey-ii+1,sizex))',(ii:min(sizey-ii+1,sizex)+ii-1)'];
        end
    end
else
    Dset = cell(sizex + sizey -1,2);
    % D0,D-1,D1,D-2,D2,D-3,D3,...D-sizex,D-(sizex+1),...
    for ii = 1:sizex
        if ii == 1
            Dset{1,1} = 0;
            Dset{1,2} = [(1:sizex)',(1:sizex)'];
        elseif ii <= sizey
            Dset{2*ii-2,1} = -(ii-1);
            Dset{2*ii-2,2} = [(ii:sizey)',(ii-1:sizey-1)'];
            Dset{2*ii - 1,1} = ii-1;
            Dset{2*ii - 1,2} = [(1:sizey-ii+1)',(ii:sizey)'];
        else
            Dset{ii + sizex - 1,1} = -(ii - 1) ;
            Dset{ii + sizex - 1,2} = [(ii:min(sizey + ii-1,sizex))',(1:min(sizey + ii-1,sizex)-ii+1)'];
        end
    end
end
%%
ksall = (-(sizex-1):sizey-1)';
k1index = (mod(ksall,3)==0);
K1 = ksall(k1index);
k2index = (mod(ksall-1,3)==0);
K2 = ksall(k2index);
k3index = (mod(ksall-2,3)==0);
K3 = ksall(k3index);

end


