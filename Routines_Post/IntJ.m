function Gs = IntJ(cCkCrd, mNDspS, mNDspE, R)
%
% Stress transformation matrix:
% T = [  c^2,  s^2,  2*s*c; ...
%        s^2,  c^2, -2*s*c; ...
%       -s*c,  s*c,  c^2-s^2 ];
%
% Strain transformation matrix:
% T = [  c^2,  s^2,  s*c; ...
%        s^2,  c^2, -s*c; ...
%       -s*c,  s*c,  c^2-s^2 ];
%
% e = [dudx, dvdy, dudy+dvdx] + e0
% s = D*e
%

global mNdCrd
global mLNodS
global cLNodE

global mGsStd_omgShp
global mGsStd_omgDrv
global vGsStd_omgWgt

global cGsEnr_omgShS
global cGsEnr_omgDvS
global cGsEnr_omgDvE
global cGsEnr_omgWgt

global constitutive_field
global constitutive_field_ddn

T2 = zeros(2,2); % coordinate axes transformation
% T3 = zeros(3,3); % stress or strain transformation

dXdx = zeros(2,2);
i0_gdv = [1,2];

nCrack = length(cCkCrd);
Gs = zeros(nCrack, 2);

for i_crk = 1:nCrack
    
    mTpCrd = cCkCrd{i_crk}([1,end],:);
    
    mTpDir(1,:) = cCkCrd{i_crk}(1,:) - cCkCrd{i_crk}(2,:);
    mTpDir(2,:) = cCkCrd{i_crk}(end,:) - cCkCrd{i_crk}(end-1,:);
    
    for i_tip = 1:2
        if R(i_crk,i_tip) > 0
            
            gs = 0;
            
            vTpCrd = mTpCrd(i_tip,:);
            vTpDir = mTpDir(i_tip,:);
            
            n = vTpDir / sqrt(vTpDir(1)^2+vTpDir(2)^2);
            
            T2(1) = n(1); T2(3) =-n(2);
            T2(2) = n(2); T2(4) = n(1);
            
            % % For stress transformation
            % T3(1) =    n(1)^2;  T3(4) =    n(2)^2;  T3(7) =   2*n(1)*n(2); 
            % T3(2) =    n(2)^2;  T3(5) =    n(1)^2;  T3(8) =  -2*n(1)*n(2);
            % T3(3) =-n(1)*n(2);  T3(6) = n(1)*n(2);  T3(9) = n(2)^2-n(1)^2;
            
            % % For strain transformation
            % T3(1) =    n(1)^2;  T3(4) =    n(2)^2;  T3(7) =     n(1)*n(2); 
            % T3(2) =    n(2)^2;  T3(5) =    n(1)^2;  T3(8) =    -n(1)*n(2);
            % T3(3) =-n(1)*n(2);  T3(6) = n(1)*n(2);  T3(9) = n(2)^2-n(1)^2;
            
            [vDmElm, mDmWgt] = Elems_sif(mNdCrd, mLNodS, vTpCrd, R(i_crk,i_tip));
            
            for i = 1:length(vDmElm)
                
                uDmElm = vDmElm(i);
                vDmWgt = mDmWgt(:,i);
                
                mElCrd = mNdCrd(mLNodS(uDmElm,:),:);
                
                if ~isempty(cLNodE{uDmElm})
                    
                    dudX = cGsEnr_omgDvS{uDmElm}*mNDspS(:,mLNodS(uDmElm,:))'*T2 ...
                         + cGsEnr_omgDvE{uDmElm}*mNDspE(:,cLNodE{uDmElm})'*T2;
                    
                    dxdX = cGsEnr_omgDvS{uDmElm}*(mElCrd-vTpCrd)*T2;
                    dqdX = cGsEnr_omgDvS{uDmElm}*vDmWgt;
                    
                    xs = cGsEnr_omgShS{uDmElm}*mElCrd;
                    qs = cGsEnr_omgShS{uDmElm}*vDmWgt;
                    
                    w = cGsEnr_omgWgt{uDmElm};
                    
                else
                    
                    dudX = mGsStd_omgDrv*mNDspS(:,mLNodS(uDmElm,:))'*T2;
                    dxdX = mGsStd_omgDrv*(mElCrd-vTpCrd)*T2;
                    dqdX = mGsStd_omgDrv*vDmWgt;
                    
                    xs = mGsStd_omgShp*mElCrd;
                    qs = mGsStd_omgShp*vDmWgt;
                    
                    w = vGsStd_omgWgt;
                    
                end
                
                i_gdv = i0_gdv;
                
                for j = 1:length(w)

                    D = constitutive_field(xs(j,:));
                    dDdx1 = constitutive_field_ddn(xs(j,:),n);
                    
                    J = dxdX(i_gdv,:);
                    detJ = det(J);
                    dA = detJ*w(j);
                    
                    dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
                    dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
                    
                    dudx = dXdx*dudX(i_gdv,:);
                    dqdx = dXdx*dqdX(i_gdv);
                    
                    e = [dudx(1); dudx(4); dudx(2)+dudx(3)];
                    
                    s = D*e;
                    
                    W = (dudx(1,1)*s(1)+dudx(1,2)*s(3))*dqdx(1) ...
                      + (dudx(1,1)*s(3)+dudx(1,2)*s(2))*dqdx(2);
                    
                    U = 0.5*(s'*e)*dqdx(1);
                    
                    gs = gs + (W - U)*dA - qs(j)*0.5*(e'*dDdx1*e)*dA;
                    % gs = gs + (W - U)*dA;
                    
                    i_gdv = i_gdv + 2;
                    
                end
                
            end
            
            Gs(i_crk,i_tip) = gs;
            
        end
    end
end
end