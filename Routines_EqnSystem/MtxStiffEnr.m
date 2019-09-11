function Kg = MtxStiffEnr(nGlDof,mNdCrd,mLNodS,vElEnr,cLNodE,mNDofE,...
                          vElPhz,cDMatx,cGsShS,cGsDvS,cGsDvE,cGsWgt)

global constitutive_field;

DIM = 2;

vElEnr = vElEnr(:)';
nElEnr = length(vElEnr);

nLNodS = size(mLNodS,2);
mLDofS = zeros(DIM,nLNodS);

nLDofS = nLNodS*DIM;
nSprGK = nElEnr*(nLDofS*5)^2; % estimate

if nSprGK == 0
    Kg = sparse([],[],[],nGlDof,nGlDof);
    return
end

vIdGdS = ones(1,nLDofS);
jSprLK = 0;

vSprGK = zeros(nSprGK,1);
jSprRw = zeros(nSprGK,1);
jSprCl = zeros(nSprGK,1);

Bs  = zeros(3,nLDofS);
iBs = DIM:DIM:nLDofS;
jBs = iBs - 1;

dXdx = zeros(2);
igd0 = 1:DIM;

for uElEnr = vElEnr
    
    % D = cDMatx{vElPhz(uElEnr)};
    
    mElCrd = mNdCrd(mLNodS(uElEnr,:),:);
    
    dNdX_s = cGsDvS{uElEnr};
    dNdX_e = cGsDvE{uElEnr};
    
    N = cGsShS{uElEnr};
    w = cGsWgt{uElEnr};
    
    vLNodE = cLNodE{uElEnr};
    nLNodE = length(vLNodE);
    nLDofE = nLNodE*DIM;
    
    Be  = zeros(3,nLDofE);
    iBe = DIM:DIM:nLDofE;
    jBe = iBe - 1;
    
    Kse = zeros(nLDofS,nLDofE);
    Kee = zeros(nLDofE);
    
    igd = igd0;
    
    for igp = 1:length(w)
        
        x = N(igp,:)*mElCrd;
        J = dNdX_s(igd,:)*mElCrd;
        
        detJ = J(1)*J(4)-J(2)*J(3);
        
        dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
        dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
        
        dNdx_s = dXdx*dNdX_s(igd,:);
        dNdx_e = dXdx*dNdX_e(igd,:);
        
        Bs(1,jBs) = dNdx_s(1,:);
        Bs(3,iBs) = dNdx_s(1,:);
        Bs(2,iBs) = dNdx_s(2,:);
        Bs(3,jBs) = dNdx_s(2,:);
        
        Be(1,jBe) = dNdx_e(1,:);
        Be(3,iBe) = dNdx_e(1,:);
        Be(2,iBe) = dNdx_e(2,:);
        Be(3,jBe) = dNdx_e(2,:);
        
        D = constitutive_field(x);
        
        DBJw = D*Be.*(detJ*w(igp));
        
        Kse = Kse + Bs'*DBJw;
        Kee = Kee + Be'*DBJw;
        
        igd = igd + DIM;
        
    end
    
    mLDofS(1,:) = DIM*mLNodS(uElEnr,:)+(1-DIM);
    for i = 2:DIM; mLDofS(i,:)=mLDofS(i-1,:)+1; end
    
    mLDofE = mNDofE(:,vLNodE);
    
    vLDofS = mLDofS(:);
    vLDofE = mLDofE(:);
    
    nElmLK = nLDofE*(2*nLDofS+nLDofE);
    jSprLK = jSprLK(end)+(1:nElmLK);
    vIdGdE = ones(1,nLDofE);
    
    mRowSE = vLDofS(:,vIdGdE); mColSE = vLDofE(:,vIdGdS)';
    mRowEE = vLDofE(:,vIdGdE); mColEE = mRowEE';
    
    jSprRw(jSprLK) = [mRowSE(:);mColSE(:);mRowEE(:)];
    jSprCl(jSprLK) = [mColSE(:);mRowSE(:);mColEE(:)];
    vSprGK(jSprLK) = [Kse(:);Kse(:);Kee(:)];

end

j = jSprLK(end)+1:nSprGK;

jSprRw(j) = [];
jSprCl(j) = [];
vSprGK(j) = [];
    
Kg = sparse(jSprRw,jSprCl,vSprGK,nGlDof,nGlDof);

end
