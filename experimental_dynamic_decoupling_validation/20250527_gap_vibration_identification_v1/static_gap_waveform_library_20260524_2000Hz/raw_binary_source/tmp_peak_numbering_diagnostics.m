clear; clc;
dataFolder = 'E:\试验数据\20260521标定数据\20260524标定_低采样率\Data';
selectedChannel = 1; bladeCount = 6; fsDefault = 2000;
detectIgnoreSeconds = 0.50; expectedBladePeriodSeconds = 2.725; eventMinDistanceFactor = 0.75; thresholdSigma = 6; minProminenceMv = 30;
files = dir(fullfile(dataFolder, 'ACTS1000_data*mm.bin'));
sourceNames = strings(numel(files),1); gaps = zeros(numel(files),1);
for k=1:numel(files)
    [~, sourceNames(k)] = fileparts(files(k).name);
    tok = regexp(char(sourceNames(k)), 'ACTS1000_data_(.+?)mm$', 'tokens', 'once');
    gaps(k)=str2double(strrep(tok{1},'_','.'));
end
[gaps,ord]=sort(gaps); files=files(ord); sourceNames=sourceNames(ord);
rows = {};
for k=1:numel(files)
    txt = fileread(fullfile(files(k).folder, sourceNames(k) + "_header.txt"));
    chanTok = regexp(txt,'chanEnableCount\s*:\s*(\d+)','tokens','once'); if isempty(chanTok), chanCount=1; else, chanCount=str2double(chanTok{1}); end
    fsTok = regexp(txt,'SampleRate\s*:\s*([-+]?\d*\.?\d+)','tokens','once'); if isempty(fsTok), fs=fsDefault; else, fs=str2double(fsTok{1}); end
    resTok = regexp(txt,'resolution\s*:\s*(\d+)','tokens','once'); if isempty(resTok), res=12; else, res=str2double(resTok{1}); end
    maxTok = regexp(txt,'0rangeMaxValue\s*:\s*([-+]?\d*\.?\d+)','tokens','once'); if isempty(maxTok), rmax=5000; else, rmax=str2double(maxTok{1}); end
    minTok = regexp(txt,'0rangeMinValue\s*:\s*([-+]?\d*\.?\d+)','tokens','once'); if isempty(minTok), rmin=-5000; else, rmin=str2double(minTok{1}); end
    fid=fopen(fullfile(files(k).folder,files(k).name),'rb'); raw=fread(fid,inf,'uint16=>uint16'); fclose(fid);
    raw=raw(1:floor(numel(raw)/chanCount)*chanCount); raw=reshape(raw,chanCount,[]);
    y=(rmax-rmin)/2^res*double(bitand(raw(selectedChannel,:),uint16(2^res-1)))+rmin;
    baseline=median(y,'omitnan'); sig=y(:)-baseline; ns=1.4826*mad(sig,1); minProm=max(thresholdSigma*ns,minProminenceMv);
    [pks,locs]=findpeaks(sig,'MinPeakDistance',round(eventMinDistanceFactor*expectedBladePeriodSeconds*fs),'MinPeakProminence',minProm);
    keep=locs>=round(detectIgnoreSeconds*fs); locs=locs(keep); pks=pks(keep);
    dt=diff(locs)/fs; if isempty(dt), medDt=NaN; cvDt=NaN; minDt=NaN; maxDt=NaN; else, medDt=median(dt); cvDt=std(dt)/mean(dt); minDt=min(dt); maxDt=max(dt); end
    rows{end+1,1}=table(sourceNames(k),gaps(k),numel(locs),medDt,cvDt,minDt,maxDt, 'VariableNames',{'sourceName','gapMm','peakCount','medianDt_s','cvDt','minDt_s','maxDt_s'});
end
T=vertcat(rows{:}); out=fullfile(dataFolder,'low_sample_gap_visualization','peak_numbering_diagnostics.csv'); writetable(T,out); disp(T); fprintf('Saved: %s\n',out);
