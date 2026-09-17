function Files=Build_R3_Fig3_NatureRelease()
%BUILD_R3_FIG3_NATURERELEASE Render the final R3 figure with SVG export.
% MATLAB exportgraphics does not accept SVG; print -dsvg preserves vector text.
root=fileparts(mfilename('fullpath')); addpath(root); addpath(fullfile(root,'local_func'),'-begin');
try
    Plot_R3_Fig3_NatureRelease;
catch ME
    if ~contains(ME.message,'svg'), rethrow(ME); end
    f=gcf; tx=findall(f,'Type','Text');
    for k=1:numel(tx)
        s=get(tx(k),'String');
        if ischar(s)&&contains(s,'Synchronous-vibration recovery'), set(tx(k),'Visible','off'); end
    end
    out=fullfile(root,'output','figures_final'); base=fullfile(out,'Fig3_R3_nature_release');
    print(f,[base '.svg'],'-dsvg'); exportgraphics(f,[base '.pdf'],'ContentType','vector'); exportgraphics(f,[base '.emf'],'ContentType','vector'); exportgraphics(f,[base '.png'],'Resolution',600);
end
Files=struct('svg',[base '.svg'],'pdf',[base '.pdf'],'emf',[base '.emf'],'png',[base '.png']); close all;
end
