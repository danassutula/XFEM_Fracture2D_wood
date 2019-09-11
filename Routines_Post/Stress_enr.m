function [s_ij,s_xy] = Stress_enr(mNDspS,mNDspE,mNdCrd,vElEnr,mLNodS,cLNodE,...
                                  vElPhz,cDMatx,mPrLod,cGsShS,cGsDvS,cGsDvE)

%--------------------------------------------------------------------------
% Stress at Gauss points (enr. el.)
%--------------------------------------------------------------------------

global constitutive_field;

vElEnr = vElEnr(:)';

ngt = 10e6; % estimate

s_ij = zeros(3,ngt);
s_xy = zeros(ngt,2);
e    = zeros(3,1);

dXdx = zeros(2,2);
igd0 = 1:2;
igt  = 0;

% for i_phz = 1:nPhase
%     
% D = cDMatx{i_phz};
%     
% for iel = vElEnr(vElPhz(vElEnr)==i_phz)
    
for iel = vElEnr
    
    mElCrd = mNdCrd(mLNodS(iel,:),:);
    
    mLDspS = mNDspS(:,mLNodS(iel,:))';
    mLDspE = mNDspE(:,cLNodE{iel})';
    
    xs   = cGsShS{iel}*mElCrd;
    dxdX = cGsDvS{iel}*mElCrd;
    dudX = cGsDvS{iel}*mLDspS+cGsDvE{iel}*mLDspE;
    
    ngp = size(xs,1);
    
    s_xy(igt+1:igt+ngp,:) = xs;
    s0 = mPrLod(:,iel);
    
    igd = igd0;
    
    for igp = 1:ngp
        
        J = dxdX(igd,:); detJ = J(1)*J(4)-J(2)*J(3);
        
        dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
        dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
        
        dudx = dXdx*dudX(igd,:);
        
        igd = igd + 2;
        igt = igt + 1;
        
        e(1) = dudx(1);
        e(2) = dudx(4);
        e(3) = dudx(2)+dudx(3);
        
        D = constitutive_field(xs(igp,:));
        
        s_ij(:,igt) = D*e + s0;
        
    end
    
end
% end

s_xy = s_xy(1:igt,:);
s_ij = s_ij(:,1:igt);

end
