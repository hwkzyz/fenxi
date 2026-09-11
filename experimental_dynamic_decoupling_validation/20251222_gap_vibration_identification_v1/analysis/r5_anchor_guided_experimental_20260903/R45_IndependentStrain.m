function T=R45_IndependentStrain()
% Independent strain frequency selection: fixed bands, no BTT frequency input.
root=fileparts(mfilename('fullpath')); pkg=fileparts(fileparts(root));
parent=fileparts(pkg); out=fullfile(root,'results','unified_closure_20260905');
alignmentFile=fullfile(parent,'20251222_low_speed_gap_prior_decoupling','outputs','Step04_BTT_STE_Resonance_Regions_20251222.mat');
A=load(alignmentFile,'strainSource'); source=A.strainSource;
tau=-source.strain_time_offset_total_sec; rows={}; nr=0;
for channel=[1 3]
 file=regexprep(source.strain_file,'AI1-03',sprintf('AI1-%02d',channel));
 S=load(file,'Datas'); t=double(S.Datas(:,1)); y=double(S.Datas(:,2));
 assert(all(diff(t)>0) && all(isfinite(y)));
 for blade=[1 5]
  if blade==1, suffix='r5_r4_surface_contract_full'; else, suffix='r5_r4_surface_b5_contract_20260904'; end
  resultFile=fullfile(pkg,'results','gap_aware',sprintf('Main_GapAware_VPTopK_FullWave_20251222_B%d_S123_%s.mat',blade,suffix));
  S=load(resultFile,'Result'); W=S.Result.WindowResult;
  for iw=1:numel(W)
   center=mean(W(iw).timeWindow)+tau; half=0.5*3*60/W(iw).rotRpmMean;
   use=t>=center-half & t<=center+half; tt=t(use)-center; yy=y(use);
   assert(numel(tt)>=20); dt=median(diff(tt));
   assert(max(abs(diff(tt)-dt))<1e-5*dt,'Strain time grid is not uniform.');
   yy=detrend(yy); N=numel(yy); nfft=2^nextpow2(N*8);
   fy=(0:floor(nfft/2))'/(dt*nfft); power=abs(fft(yy.*hann(N),nfft)); power=power(1:numel(fy));
   bands=[300 1000;520 640];
   for ib=1:2
    band=bands(ib,:); candidates=find(fy>=band(1)&fy<=band(2));
    [~,ip]=max(power(candidates)); f0=fy(candidates(ip));
    lo=max(band(1),f0-2/(N*dt)); hi=min(band(2),f0+2/(N*dt));
    % Coarse scan avoids assuming the padded FFT peak is a unique LS minimum.
    grid=lo:0.1:hi; costs=arrayfun(@(f)tone_cost(f,tt,yy),grid); [~,k]=min(costs);
    low=grid(max(1,k-1)); high=grid(min(numel(grid),k+1));
    f=fminbnd(@(f)tone_cost(f,tt,yy),low,high,optimset('TolX',1e-5));
    [sse,amp]=tone_cost(f,tt,yy); fit=W(iw).modelFits.gap_only;
    nr=nr+1; rows{nr}=table(blade,channel,iw,band(1),band(2),center,2*half,N, ...
     f,amp,1-sse/sum(yy.^2),1/(N*dt),fit.freqHz,fit.EO,W(iw).rotRpmMean/60, ...
     fit.freqHz-f,round(f/(W(iw).rotRpmMean/60)), ...
     'VariableNames',{'blade','channel','window','band_low_hz','band_high_hz', ...
     'strain_center_sec','duration_sec','samples','strain_frequency_hz','strain_amplitude_raw', ...
     'single_tone_explained_fraction','rayleigh_resolution_hz','btt_frequency_hz','btt_eo', ...
     'rot_frequency_hz','btt_minus_strain_hz','strain_nearest_order'});
   end
  end
 end
end
T=vertcat(rows{:}); writetable(T,fullfile(out,'IndependentStrain_PerWindow.csv'));
save(fullfile(out,'IndependentStrain.mat'),'T','alignmentFile','source','tau');
Summary=groupsummary(T,{'blade','channel','band_low_hz'},{'mean','std','min','max'}, ...
 {'btt_minus_strain_hz','single_tone_explained_fraction'});
writetable(Summary,fullfile(out,'IndependentStrain_Summary.csv')); disp(Summary);
end
function [sse,amp]=tone_cost(f,t,y)
H=[cos(2*pi*f*t),sin(2*pi*f*t),ones(size(t))]; b=H\y; e=y-H*b;
sse=e'*e; amp=hypot(b(1),b(2));
end
