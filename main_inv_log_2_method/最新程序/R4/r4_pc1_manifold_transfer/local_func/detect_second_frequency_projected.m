function info = detect_second_frequency_projected(map, singleFit, lib, cfg)
%DETECT_SECOND_FREQUENCY_PROJECTED Noise-aware residual-order diagnostic.
% The complete trusted waveform is used. Candidate second-frequency columns
% are projected off gap, translation, and the fitted first-frequency tangent
% space before their conditional reduction in voltage SSE is measured.

t=map.t_v(:);x=map.x_v(:);V=map.V_a(:);t0=mean(t);
g=singleFit.g;dx=singleFit.dx;A=singleFit.A;phi=singleFit.phi;f1=singleFit.f;
s1=sin(2*pi*f1*t+phi);c1=cos(2*pi*f1*t+phi);
z=x-dx-A*s1;
V1=eval_gap_template(lib,g,z);Fx=eval_gap_derivative(lib,g,z);
h=get_field(cfg,'informationGapStepMm',1e-3);
Fg=(eval_gap_template(lib,g+h,z)-eval_gap_template(lib,g-h,z))/(2*h);
r=V-V1;
Jeta=[Fg,-Fx,-Fx.*s1,-Fx.*c1,-Fx.*A.*c1.*(2*pi*(t-t0))];
valid=isfinite(r)&all(isfinite(Jeta),2)&isfinite(Fx);
t=t(valid);r=r(valid);Fx=Fx(valid);Jeta=Jeta(valid,:);
n=numel(r);
info=empty_info(f1,n);
if n<=size(Jeta,2)+4,return;end

scale=sqrt(sum(Jeta.^2,1));keep=scale>max(scale)*1e-10;
Jeta=Jeta(:,keep)./max(scale(keep),eps);
[Q,R]=qr(Jeta,0);
rk=nnz(abs(diag(R))>max(abs(diag(R)))*1e-10);
Q=Q(:,1:rk);r=r-Q*(Q.'*r);
sse1=dot(r,r);sigma2=sse1/max(n-rk,1);

grid=get_field(cfg,'route30SingleGrid',get_field(cfg,'f1Grid',300:5:1500));
grid=unique(grid(:).');
minSep=get_field(cfg,'route30SecondMinSeparationHz',20);
grid=grid(abs(grid-f1)>=minSep);
if isempty(grid)||sse1<=eps,return;end

S=sin(2*pi*t*grid);C=cos(2*pi*t*grid);
QS=-Fx.*S;QC=-Fx.*C;
if rk>0
    QS=QS-Q*(Q.'*QS);QC=QC-Q*(Q.'*QC);
end
gss=sum(QS.^2,1);gcc=sum(QC.^2,1);gsc=sum(QS.*QC,1);
bs=r.'*QS;bc=r.'*QC;detG=gss.*gcc-gsc.^2;
usable=detG>max(detG)*1e-12;
gain=zeros(size(grid));a=zeros(size(grid));b=zeros(size(grid));
a(usable)=(gcc(usable).*bs(usable)-gsc(usable).*bc(usable))./detG(usable);
b(usable)=(gss(usable).*bc(usable)-gsc(usable).*bs(usable))./detG(usable);
gain(usable)=bs(usable).*a(usable)+bc(usable).*b(usable);
gain=max(gain,0);[bestGain,k]=max(gain);

sse2=max(sse1-bestGain,eps);
deltaK=3;
deltaBic=n*log(sse2/sse1)+deltaK*log(n);
amp=hypot(a(k),b(k));
if amp>0&&usable(k)
    Cbeta=sigma2/detG(k)*[gcc(k),-gsc(k);-gsc(k),gss(k)];
    direction=[a(k);b(k)]/amp;
    ampStd=sqrt(max(direction.'*Cbeta*direction,0));
else
    ampStd=Inf;
end
ampZ=amp/max(ampStd,eps);
strong=deltaBic<get_field(cfg,'route30SecondDeltaBicThreshold',-10)&&...
    ampZ>get_field(cfg,'route30SecondAmplitudeZThreshold',3);
info=struct('tested',true,'first_frequency_hz',f1,...
    'second_frequency_hz',grid(k),'conditional_amplitude_mm',amp,...
    'conditional_amplitude_std_mm',ampStd,'amplitude_z',ampZ,...
    'projected_sse_single',sse1,'projected_sse_dual',sse2,...
    'sse_reduction',bestGain,'delta_bic',deltaBic,...
    'n_samples',n,'nuisance_rank',rk,'strong_second_frequency',strong);
end

function info=empty_info(f1,n)
info=struct('tested',false,'first_frequency_hz',f1,...
    'second_frequency_hz',NaN,'conditional_amplitude_mm',NaN,...
    'conditional_amplitude_std_mm',Inf,'amplitude_z',0,...
    'projected_sse_single',NaN,'projected_sse_dual',NaN,...
    'sse_reduction',0,'delta_bic',Inf,'n_samples',n,...
    'nuisance_rank',0,'strong_second_frequency',false);
end

function v=get_field(s,name,d)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)),v=s.(name);else,v=d;end
end
