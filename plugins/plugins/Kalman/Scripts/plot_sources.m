full_address = 'exportImage';
file_name = 'carsten';
file_full = fullfile(full_address, file_name);

%%

%% Deep Component
set(zef.h_axes1,'clim',[0 1])
f1 = gcf;
f2 = figure();
%612   510
f2.Position(3) = 620;
f2.Position(4) = 600;
copyobj(findobj(f1.Children, 'tag','axes1'), f2);

f2.Children.Position(2) = 60;
grid off
axis off
%view(-90,40)
view(-120,25)
f2.Children.Children(end).FaceAlpha = 0.45;
f2.Children.Children(end-1).FaceAlpha = 0.45;


set(gcf, 'color', 'w')
%%
print(f2,'-dpng','-r150',[file_full 'frontDeep' '.png'])
savefig(f2, [file_full 'frontDeep' '.fig'])


view(-85,5)
print(f2,'-dpng','-r150',[file_full 'sideDeep' '.png'])
savefig(f2, [file_full 'sideDeep' '.fig'])

%%

%% Cortex component
set(zef.h_axes1,'clim',[0 1])
f1 = gcf;
f2 = figure();
copyobj(findobj(f1.Children, 'tag','axes1'), f2);

f2.Children.Position(2) = 35;
grid off
axis off
view(-140,60)

%f2.Children.View = [-110 40];

set(gcf, 'color', 'w')
f2.Children.Children(end).FaceAlpha = 0.35;
f2.Children.Children(end-1).FaceAlpha = 0.35;
f2.Children.Children(end-2).FaceAlpha = 0.25;


%%
view(0,0)
print(f2,'-dpng','-r150',[file_full 'cortex1' '.png'])
savefig(f2, [file_full 'cortex1' '.fig'])

view(-90,0)
print(f2,'-dpng','-r150',[file_full 'cortex2' '.png'])
savefig(f2, [file_full 'cortex2' '.fig'])

view(90,0)
print(f2,'-dpng','-r150',[file_full 'cortex3' '.png'])
savefig(f2, [file_full 'cortex3' '.fig'])

