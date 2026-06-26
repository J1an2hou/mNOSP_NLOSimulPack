function [H_f,dH_fdk,d2H_fdk2] ...
    = Hk_TimeAve_v2(nq,NT,nw,nR,omega,A0,eta,R,HR,wtR,wc,kpt,alatt,blatt)
Hk_timeq = zeros(nw,nw,nq);
Hk_time_q = zeros(nw,nw,nq);
dHdk_timeq = zeros(nw,nw,3,nq);
d2Hdk2_timeq = zeros(nw,nw,6,nq);
dHdk_time_q = zeros(nw,nw,3,nq);
d2Hdk2_time_q = zeros(nw,nw,6,nq);

hr0 = FloPhase_ave(omega,A0,eta,R,wc,alatt,nw,nR,0,NT);
[Hk_q0, dHdk_q0, d2Hdk2_q0] = buildHk(R,HR.*hr0,wtR,wc,nw,nR,kpt,alatt,blatt);
for iq = 1:nq
    hrq = FloPhase_ave(omega,A0,eta,R,wc,alatt,nw,nR,iq,NT);
    hr_q = FloPhase_ave(omega,A0,eta,R,wc,alatt,nw,nR,-iq,NT);
    [Hk_timeq(:,:,iq), dHdk_timeq(:,:,:,iq), d2Hdk2_timeq(:,:,:,iq)] ...
        = buildHk(R,HR.*hrq,wtR,wc,nw,nR,kpt,alatt,blatt);
    [Hk_time_q(:,:,iq), dHdk_time_q(:,:,:,iq), d2Hdk2_time_q(:,:,:,iq)] ...
        = buildHk(R,HR.*hr_q,wtR,wc,nw,nR,kpt,alatt,blatt);
end

H_f = Hk_q0;
dH_fdk = dHdk_q0;
d2H_fdk2 = d2Hdk2_q0;
for iq = 1:nq
    H_f = H_f + comm(Hk_time_q(:,:,iq),Hk_timeq(:,:,iq),-1)/omega/iq;
    for i = 1:3
        dH_fdk(:,:,i) = dH_fdk(:,:,i) ...
            + comm(dHdk_time_q(:,:,i,iq),Hk_timeq(:,:,iq),-1)/omega/iq ...
            + comm(Hk_time_q(:,:,iq),dHdk_timeq(:,:,i,iq),-1)/omega/iq;
    end
    for j = 1:6
        [ja,jb] = invVoigt(j);
        d2H_fdk2(:,:,j) = d2H_fdk2(:,:,j) ...
            + comm(d2Hdk2_time_q(:,:,j,iq),Hk_timeq(:,:,iq),-1)/omega/iq ...
            + comm(dHdk_time_q(:,:,ja,iq),dHdk_timeq(:,:,jb,iq),-1)/omega/iq ...
            + comm(dHdk_time_q(:,:,jb,iq),dHdk_timeq(:,:,ja,iq),-1)/omega/iq ...
            + comm(Hk_time_q(:,:,iq),d2Hdk2_timeq(:,:,j,iq),-1)/omega/iq;
    end
end
H_f = roundn(H_f, -6);
if ~ishermitian(H_f)
    H_f = Hert_mat(H_f);
end
for i = 1:3
    dH_fdk(:,:,i) = roundn(dH_fdk(:,:,i), -6);
    if ~ishermitian(dH_fdk(:,:,i))
        dH_fdk(:,:,i) = Hert_mat(dH_fdk(:,:,i));
    end
end
for j = 1:6
    d2H_fdk2(:,:,j) = roundn(d2H_fdk2(:,:,j), -6);
    if ~ishermitian(d2H_fdk2(:,:,j))
        d2H_fdk2(:,:,j) = Hert_mat(d2H_fdk2(:,:,j));
    end
end

end