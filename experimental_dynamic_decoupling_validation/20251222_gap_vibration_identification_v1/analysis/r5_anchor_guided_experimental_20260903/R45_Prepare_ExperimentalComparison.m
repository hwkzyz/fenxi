function Manifest=R45_Prepare_ExperimentalComparison()
% Frozen B2 registration, target-excluded PCA, common-sensor profile objective.
root=fileparts(mfilename('fullpath')); pkg=fileparts(fileparts(root));
out=fullfile(root,'results','unified_closure_20260905');
S=load(fullfile(root,'results','r4_frontend_only','r4_experimental_surface_family.mat'),'F'); F=S.F;
S=load(fullfile(root,'results','r4_frontend_only','r4_b2_anchor_registration.mat'),'G'); G=S.G;
S=load(fullfile(pkg,'inputs','prepared','foundation','step04_low_speed_template', ...
 'Template_OPRCenterStd_AdaptiveSG_LowSpeed_AllBlades_S123_20251222.mat'),'Template'); Template=S.Template;
rows={}; nr=0;
for target=[1 5]
 keep=F.blade_id~=target; Q=F.q(:,keep); mu=mean(Q,2); [U,~,~]=svd(Q-mu,'econ');
 [z,ord]=sort(U(:,1)'*(Q-mu)); Q=Q(:,ord); zg=linspace(z(1),z(end),301);
 C=interp1(z,Q.',zg,'linear').'; gg=linspace(min(F.gap_mm),max(F.gap_mm),281);
 Xg=[ones(1,numel(gg));1./gg;log(gg/F.gref_mm)];
 total=zeros(301,1); count=0; independent=[]; profiles=[];
 for sid=[1 2 3]
  a=Template.SensorBlade([Template.SensorBlade.blade_id]==target & [Template.SensorBlade.sensor_id]==sid);
  r=G.registration([G.registration.sensor_id]==sid); xa=double(a.x_grid(:)); ya=(double(a.v_grid(:))-double(a.baseline))*1000;
  m=isfinite(xa)&isfinite(ya)&logical(a.valid_grid_mask(:))&logical(a.domain_effective_mask(:));
  y=interp1(r.x_scale*(xa(m)-r.tau_mm),ya(m),F.x_mm,'pchip',NaN);
  use=F.support_mask(:)&isfinite(y); yy=y(use)-r.voltage_offset_mv;
  J=zeros(301,numel(gg));
  for iz=1:301
   B=reshape(C(:,iz),numel(F.x_mm),3); B=r.voltage_gain*B(use,:);
   J(iz,:)=sqrt(max(0,sum(Xg.*((B'*B)*Xg),1)-2*(yy'*B)*Xg+yy'*yy)/nnz(use));
  end
  profile=min(J,[],2); [~,ii]=min(profile); independent(end+1)=zg(ii); %#ok<AGROW>
  profiles(:,sid)=profile; total=total+nnz(use)*profile.^2; count=count+nnz(use); %#ok<AGROW>
 end
 objective=sqrt(total/count); [best,iz]=min(objective);
 spread=range(independent)/range(z); gate=spread<=0.10;
 S=load(fullfile(pkg,'inputs','calibration',sprintf('Step06I_OffsetTiltShared_GapLibrary_20251222_B%d_S123.mat',target)),'CorrectedGapLibrary');
 base=S.CorrectedGapLibrary;
 methods={'Oracle','CommonB2','R5LOBO'};
 columns={F.q(:,F.blade_id==target),F.q(:,F.blade_id==2),C(:,iz)};
 for k=1:3
  CorrectedGapLibrary=base; CorrectedGapLibrary.responseSurface.coeff=reshape(columns{k},numel(F.x_mm),3);
  CorrectedGapLibrary.responseSurface.method=['R45 closure ' methods{k}];
  CorrectedGapLibrary.r45=struct('target',target,'method',methods{k},'training_blades',F.blade_id(keep), ...
   'common_latent',zg(iz),'sensor_latents',independent,'sensor_spread_fraction',spread,'gate_pass',gate, ...
   'anchor_rmse_mv',best,'diagnostic_if_gate_fails',true);
  assert(isequaln(CorrectedGapLibrary.sensor,base.sensor));
  rb=base.responseSurface; rc=CorrectedGapLibrary.responseSurface;
  assert(isequaln(rmfield(rb,{'coeff','method'}),rmfield(rc,{'coeff','method'})));
  file=fullfile(out,sprintf('R45_%s_B%d.mat',methods{k},target)); save(file,'CorrectedGapLibrary','-v7.3');
  nr=nr+1; rows{nr}=table(target,string(methods{k}),string(file),spread,gate,best, ...
   'VariableNames',{'target_blade','method','library_file','sensor_spread_fraction','gate_pass','anchor_rmse_mv'});
 end
 save(fullfile(out,sprintf('LOBO_Localization_B%d.mat',target)),'zg','gg','profiles','objective','independent','keep');
end
Manifest=vertcat(rows{:}); writetable(Manifest,fullfile(out,'ExperimentalComparison_Manifest.csv')); disp(Manifest);
end
