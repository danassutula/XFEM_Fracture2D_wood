function Kg = MtxStiffStd(nGlDof,mNdCrd,mElNod,vElPhz,cDMatx,dNdX,N,w)

global constitutive_field;

DIM = 2;

[nElems,nElNod] = size(mElNod);
mElDof = zeros(DIM,nElNod);

nElDof = nElNod*DIM;
nElmLK = nElDof^2;
nSprGK = nElmLK*nElems;

if nSprGK == 0
    Kg = sparse([],[],[],nGlDof,nGlDof);
    return
end

Ke0 = zeros(nElDof,nElDof);

jSprLK = 1:nElmLK;
vIdGrd = ones(1,nElDof);

vSprGK = zeros(nSprGK,1);
jSprRw = zeros(nSprGK,1);
jSprCl = zeros(nSprGK,1);

B = zeros(3,nElDof);
iB = DIM:DIM:nElDof;
jB = iB - 1;

igd0 = 1:DIM;
igp0 = 1:length(w);
dXdx = zeros(2,2);

for iElemn = 1:nElems
    
    % D = cDMatx{vElPhz(iElemn)};
    
    mElCrd = mNdCrd(mElNod(iElemn,:),:);
    
    Ke = Ke0;
    igd = igd0;
    
    for igp = igp0

        x = N(igp,:)*mElCrd;
        J = dNdX(igd,:)*mElCrd;
        
        detJ = J(1)*J(4)-J(2)*J(3);

        dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
        dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;

        dNdx = dXdx*dNdX(igd,:);

        B(1,jB) = dNdx(1,:);
        B(3,iB) = dNdx(1,:);
        B(2,iB) = dNdx(2,:);
        B(3,jB) = dNdx(2,:);

        D = constitutive_field(x);
        
        Ke = Ke + B'*D*B.*(detJ*w(igp));

        igd = igd + DIM;

    end

    mElDof(1,:) = mElNod(iElemn,:) .* DIM + (1-DIM);
    for i = 2:DIM; mElDof(i,:) = mElDof(i-1,:) + 1; end

    vElDof = mElDof(:);
    mRwGrd = vElDof(:,vIdGrd);
    mClGrd = mRwGrd';

    jSprRw(jSprLK) = mRwGrd(:);
    jSprCl(jSprLK) = mClGrd(:);
    vSprGK(jSprLK) = Ke(:);

    jSprLK = jSprLK + nElmLK;

end

Kg = sparse(jSprRw,jSprCl,vSprGK,nGlDof,nGlDof);

end
