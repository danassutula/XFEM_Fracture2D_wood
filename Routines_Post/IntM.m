function Gs_mix = IntM(cCkCrd, mNDspS, mNDspE, R, mu_tip, kappa_tip, K1_aux, K2_aux)
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

s = zeros(2,2); % std. state
S = zeros(2,2); % aux. state

e = zeros(3,1);
dXdx = zeros(2,2);
i0_gdv = [1,2];

nCrack = length(cCkCrd);
Gs_mix = zeros(nCrack, 2);

for i_crk = 1:nCrack
    
    mTpCrd = cCkCrd{i_crk}([1,end],:);
    
    mTpDir(1,:) = cCkCrd{i_crk}(1,:)   - cCkCrd{i_crk}(2,:);
    mTpDir(2,:) = cCkCrd{i_crk}(end,:) - cCkCrd{i_crk}(end-1,:);
    
    for i_tip = 1:2
        if R(i_crk,i_tip) > 0
            
            gs_mix = 0;
            % gs_st1 = 0;
            % gs_st2 = 0;
            
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
                mElCrd_local = (mElCrd-vTpCrd)*T2;
                
                if ~isempty(cLNodE{uDmElm})
                    
                    dudX = cGsEnr_omgDvS{uDmElm}*mNDspS(:,mLNodS(uDmElm,:))'*T2 ...
                         + cGsEnr_omgDvE{uDmElm}*mNDspE(:,cLNodE{uDmElm})'*T2;
                    
                    dxdX = cGsEnr_omgDvS{uDmElm}*mElCrd_local;
                    dqdX = cGsEnr_omgDvS{uDmElm}*vDmWgt;
                    
                    xs = cGsEnr_omgShS{uDmElm}*mElCrd_local;
                    qs = cGsEnr_omgShS{uDmElm}*vDmWgt;
                    
                    xs_global = cGsEnr_omgShS{uDmElm}*mElCrd;
                    
                    w = cGsEnr_omgWgt{uDmElm};
                    
                else
                    
                    dudX = mGsStd_omgDrv*mNDspS(:,mLNodS(uDmElm,:))'*T2;
                    dxdX = mGsStd_omgDrv*mElCrd_local;
                    dqdX = mGsStd_omgDrv*vDmWgt;
                    
                    xs = mGsStd_omgShp*mElCrd_local;
                    qs = mGsStd_omgShp*vDmWgt;
                    
                    xs_global = mGsStd_omgShp*mElCrd;
                
                    w = vGsStd_omgWgt;
                    
                end
                
                %----------------------------------------------------------
                % Auxiliary state
                %----------------------------------------------------------
                
                [theta,radii] = cart2pol(xs(:,1), xs(:,2));
                
                [~,DUDX,DUDY] = CrackTipField_Disp(...
                    K1_aux, K2_aux, mu_tip, kappa_tip, radii, theta);
                
                [SXX,SYY,SXY] = CrackTipField_Stress(...
                    K1_aux, K2_aux, radii, theta);
                
                EXX = DUDX(:,1);
                EYY = DUDY(:,2);
                EXY = DUDX(:,2) + DUDY(:,1);
                
                %----------------------------------------------------------
                
                i_gdv = i0_gdv;
                
                for j = 1:length(w)
                    
                    D = constitutive_field(xs_global(j,:));
                    dDdx1 = constitutive_field_ddn(xs_global(j,:),n);
                    
                    J = dxdX(i_gdv,:);
                    detJ = det(J);
                    dA = detJ*w(j);
                    
                    dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
                    dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
                    
                    dudx = dXdx*dudX(i_gdv,:);
                    dqdx = dXdx*dqdX(i_gdv);
                    
                    exx = dudx(1);
                    eyy = dudx(4);
                    exy = dudx(2)+dudx(3);
                    
                    sxx = D(1)*exx + D(4)*eyy + D(7)*exy;
                    syy = D(2)*exx + D(5)*eyy + D(8)*exy;
                    sxy = D(3)*exx + D(6)*eyy + D(9)*exy;
                    
                    s(1) = sxx; s(3) = sxy;
                    s(2) = sxy; s(4) = syy;
                    
                    S(1) = SXX(j); S(3) = SXY(j);
                    S(2) = SXY(j); S(4) = SYY(j);
                    
                    %------------------------------------------------------
                    %  Mixed State
                    %------------------------------------------------------
                    
                    W = (dudx(1,:)*S+DUDX(j,:)*s)*dqdx;
                    
                    U = 0.5*(exx*SXX(j) + eyy*SYY(j) + exy*SXY(j) + ...
                             EXX(j)*sxx + EYY(j)*syy + EXY(j)*sxy)*dqdx(1);
                    
                    U_D = qs(j)*([exx,eyy,exy]*dDdx1*[EXX(j);EYY(j);EXY(j)]);
                    
                    gs_mix = gs_mix + (W - U - U_D)*dA;
                    % gs_mix = gs_mix + (W - U)*dA; 
                    
                    i_gdv = i_gdv + 2;
                    
                end
                
            end
            
            Gs_mix(i_crk,i_tip) = gs_mix;
            
        end
    end
end
end