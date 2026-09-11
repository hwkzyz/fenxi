function T = R45_Validate_UnifiedTransfer()
% Same no-extrapolation PC1 coefficient interpolation and gap profiling in
% simulation and measured-static leave-surface-out tests. All errors in mV.
root=fileparts(mfilename('fullpath')); pkg=fileparts(fileparts(root));
workspace=fileparts(fileparts(pkg));
out=fullfile(root,'results','unified_closure_20260905');
if ~isfolder(out), mkdir(out); end
S=load(fullfile(workspace,'main_inv_log_2_method','最新程序', ...
 '33_tilt_calibration_transfer','output','waveform_families','TiltGap_WaveformSurfaceData.mat'),'Data');
D=S.Data;
datasets={struct('name','R4','Y',permute(double(D.response_V)*1000,[1 3 2]), ...
 'g',D.gap_mm(:),'gref',0.8,'support',true(numel(D.x_mm),1)),[]};
S=load(fullfile(root,'results','r4_frontend_only','r4_experimental_surface_family.mat'),'F'); F=S.F;
datasets{2}=struct('name','R5_static','Y',F.waveforms_mv,'g',F.gap_mm(:), ...
 'gref',F.gref_mm,'support',F.support_mask(:));
rows={}; n=0;
for id=1:numel(datasets)
 D=datasets{id}; Y=D.Y; use=D.support; g=D.g; nx=size(Y,1); nb=size(Y,2); ng=numel(g);
 assert(all(isfinite(Y(:))) && all(g>0));
 X=[ones(ng,1),1./g,log(g/D.gref)];
 for target=1:nb
  train=setdiff(1:nb,target); Q=zeros(3*nx,numel(train));
  for k=1:numel(train), B=(X\squeeze(Y(:,train(k),:)).').'; Q(:,k)=B(:); end
  mu=mean(Q,2); [U,~,~]=svd(Q-mu,'econ'); z=U(:,1)'*(Q-mu); [z,ord]=sort(z); Q=Q(:,ord);
  zg=linspace(z(1),z(end),301); Qi=interp1(z,Q.',zg,'linear').';
  % Methods share training folds, support, gap search, and anchor objective.
  candidates={mu,Q,Qi}; names={'Common','Nearest','PC1'};
  for ia=1:ng
   y=Y(use,target,ia); other=setdiff(1:ng,ia); truth=squeeze(Y(use,target,:));
   for unknown=[false true]
    if unknown, gs=linspace(min(g),max(g),281); else, gs=g(ia); end
    Xs=[ones(1,numel(gs));1./gs;log(gs/D.gref)];
    for method=1:3
     C=candidates{method}; J=zeros(size(C,2),numel(gs));
     for iz=1:size(C,2)
      B=reshape(C(:,iz),nx,3); B=B(use,:);
      % Algebraically exact SSE via the 3-column Gram matrix.
      J(iz,:)=sqrt(max(0,sum(Xs.*((B'*B)*Xs),1)-2*(y'*B)*Xs+y'*y)/nnz(use));
      if iz==1
       direct=sqrt(mean((B*Xs-y).^2,1));
       assert(max(abs(direct-J(iz,:)))<1e-4,'Gram SSE differs from direct RMSE.');
      end
     end
     [profile,igs]=min(J,[],2); [anchorErr,iz]=min(profile);
     B=reshape(C(:,iz),nx,3); pred=B(use,:)*X.'; e=pred(:,other)-truth(:,other);
     surfaceErr=sqrt(mean(e(:).^2)); normRange=range(truth(:,other),'all');
     % Measured gap increments relative to each anchor, evaluated at held gaps.
     incErr=(pred(:,other)-pred(:,ia))-(truth(:,other)-truth(:,ia));
     near=find(profile<=anchorErr+1); worst=0; spread=0;
     for j=near(:).'
      Bj=reshape(C(:,j),nx,3); P=Bj(use,:)*X.'; E=P(:,other)-truth(:,other);
      worst=max(worst,sqrt(mean(E(:).^2)));
      E=P(:,other)-pred(:,other); spread=max(spread,sqrt(mean(E(:).^2)));
     end
     n=n+1; rows{n}=table(string(D.name),target,ia,unknown,string(names{method}), ...
      anchorErr,surfaceErr,100*surfaceErr/normRange,sqrt(mean(incErr(:).^2)), ...
      numel(near),worst,spread,gs(igs(iz)),g(ia),iz==1||iz==size(C,2), ...
      'VariableNames',{'dataset','target_index','anchor_index','hidden_gap','method', ...
      'anchor_rmse_mv','held_surface_rmse_mv','held_surface_nrmse_pct', ...
      'gap_increment_rmse_mv','near_count','near_worst_surface_rmse_mv', ...
      'near_max_distance_from_selected_mv','fitted_gap_mm','true_anchor_gap_mm','candidate_boundary'});
    end
   end
  end
  fprintf('%s target %d/%d done\n',D.name,target,nb);
 end
end
T=vertcat(rows{:}); writetable(T,fullfile(out,'UnifiedTransfer_PerAnchor.csv'));
Summary=groupsummary(T,{'dataset','hidden_gap','method'}, {'median','max'}, ...
 {'held_surface_nrmse_pct','gap_increment_rmse_mv','near_max_distance_from_selected_mv'});
writetable(Summary,fullfile(out,'UnifiedTransfer_Summary.csv'));
save(fullfile(out,'UnifiedTransfer.mat'),'T','Summary','-v7.3'); disp(Summary);
end
