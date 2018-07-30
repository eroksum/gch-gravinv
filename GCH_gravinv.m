function GCH_gravinv
%% Gravity anomaly inversion based on Granser–Cordell and Henderson’s 
% algorithm
%%Code by: Luan Thanh Pham,Erdinc Oksum and Thanh Doc Do
%%email: eroksum@gmail.com
clc;clear all;clear global;
%%%%%%%%%%%%%%%%%%%%%%%%%%%Structure of main figure window
fig1 = figure('MenuBar','none','Name',...
'GCHinv','NumberTitle','off','ReSize','off',...
'Tag','mainfig','Color','w',...
'units','normalized','outerposition',[0 0 1 1]);
%%%%Menu Items 
menu1 = uimenu('Parent',fig1,'Label','IMPORT DATA','Tag','men1'); 
uimenu('Parent',menu1,'Label','2-D Grid (*.grd)','CallBack',@impgrd);
uimenu('Parent',menu1,'Label','XYZ grid (*.dat)','CallBack',@impgrdat);
menu3=uimenu('Parent',fig1,'Label','SAVE','Tag','savbut','enable','off');
uimenu('Parent',menu3,'Label','Export ascii grid (*.grd)','CallBack',@sav1data);
uimenu('Parent',menu3,'Label','Export mat file (*.mat)','CallBack',@sav2data);
uimenu('Parent',fig1,'Tag','expbut','enable','off','Label','Export_Fig.','CallBack',@expfig);
uimenu('Parent',fig1,'Label','Forward Calc','CallBack',@forwardc);
%%%% Items in Window
s1='GCH_gravinv : Gravity anomaly inversion based on Granser–Cordell and Henderson’s algorithm';
s2='Code by : Luan Thanh Pham, Erdinc Oksum and Thanh Duc Do';
uicontrol('Parent',fig1,'Style','text','units','normalized','Position',...
[0.2 0.01 0.8 0.03],'FontWeight','normal','HorizontalAlignment','left',...
'BackGroundColor','w','ForeGroundColor',[.2 .2 .2],...
'string',[s1 '--' s2])
%%%%panels on lefthand side
panel1=uipanel('Parent',fig1,'Position',[0,0,0.2,1],...
'BackgroundColor',[.90 .90 .90],'Tag','pan1','BorderType','line',...
'Title','Control Panel','ForegroundColor','k','TitlePosition','lefttop');
panel2=uipanel('Parent',panel1,'Position',[0.05,.6,0.9,.4],...
'BackgroundColor',[.96 .96 .96],'Tag','pan2','BorderType','beveledin',...
'Title','Grid Info & Settings','ForegroundColor','k','TitlePosition','centertop');
uicontrol('Parent',panel2,'Style','listbox','units','normalized','Position',[.02,.7,.96,.28],...
'String',{''},'Tag','listb1')
%Actual grid spacing in km
uicontrol('Parent',panel2,'Style','text','units','normalized',...
'Position',[0.02 0.61 0.65 0.08],...
'FontWeight','bold','HorizontalAlignment','left','BackGroundColor',[.96 .96 .96],...
'ForeGroundColor','k','string','Actual grid spacing (in km):')
uicontrol('Parent',panel2,'Style','edit','units','normalized',...
'Position',[0.67 0.63 0.31 0.06],...
'FontWeight','bold','HorizontalAlignment','center','BackGroundColor',[.90 .90 .90],...
'ForeGroundColor','r','string','1','Tag','edactdx')
%r0: the density contrast observed at the ground surface
uicontrol('Parent',panel2,'Style','text','units','normalized',...
'Position',[0.02 0.54 0.55 0.06],...
'FontWeight','bold','HorizontalAlignment','left','BackGroundColor',[.96 .96 .96],...
'ForeGroundColor','k','string','density contrast r0:',...
'TooltipString','the density contrast observed at the ground surface')
uicontrol('Parent',panel2,'Style','edit','units','normalized',...
'Position',[0.6 0.54 0.38 0.06],...
'FontWeight','bold','HorizontalAlignment','center','BackGroundColor',[.90 .90 .90],...
'ForeGroundColor','k','string','-0.0','Tag','edro')
%lambda: the decay constant
uicontrol('Parent',panel2,'Style','text','units','normalized',...
'String','decay constant (lambda)','Position',[0.02 .43 .55 0.1],...
'TooltipString','density contrast varies exponentially with depth, define the decay constant',...
'FontWeight','bold','BackGroundColor',[.96 .96 .96],...
'HorizontalAlignment','left')
uicontrol('Parent',panel2,'Style','edit','units','normalized',...
'Position',[0.6 0.45 0.38 0.06],...
'FontWeight','bold','HorizontalAlignment','center','BackGroundColor',[.90 .90 .90],...
'ForeGroundColor','k','string','0','Tag','edlambda')
%%%%%%ignoring data from edges (for define the area for rms calculating)
uicontrol('Parent',panel2,'Style','text','units','normalized',...
'String','Ignoring from edges','Position',[0.02 .32 .55 0.06],...
'TooltipString','n nodes to be ignored from edges of the data matrix, default is 0 for all area',...
'FontWeight','bold','BackGroundColor',[.96 .96 .96],...
'HorizontalAlignment','left')
uicontrol('Parent',panel2,'Style','edit','units','normalized',...
'Position',[0.6 0.32 0.38 0.06],...
'FontWeight','bold','HorizontalAlignment','center','BackGroundColor',[.90 .90 .90],...
'ForeGroundColor','k','string','0','Tag','edign')
%%iteration criteria
bg1=uibuttongroup(panel2,'units','normalized','Position',[0.01 0.005 0.979 0.3],...
'Title','Iteration Stop','TitlePosition','lefttop','SelectionChangeFcn',@radbehav);
%%%radiobuttons in button group
uicontrol(bg1,'Style','radiobutton','String','divergence [ RMS(i+1) > RMS(i) ]',...
'units','normalized','Position',[.01 .65 .9 .3],'Tag','rbut1');
uicontrol(bg1,'Style','radiobutton','String','convergence [ RMS < threshold ]',...
'units','normalized', 'Position',[.01 .25 .9 .3],'Tag','rbut2');
%%%%threshold edit box
uicontrol('Parent',bg1,'Style','edit','units','normalized',...
'Position',[0.78 0.25 0.20 0.3],...
'FontWeight','bold','HorizontalAlignment','center','BackGroundColor',[.90 .90 .90],...
'ForeGroundColor','k','string','0.01','Tag','edcriterio','enable','off')
%%%start button of inversion
uicontrol('Parent',panel1,'Style','pushbutton','units','normalized',...
'Position',[0.5 0.54 0.45 0.05],'ForeGroundColor','b',...
'FontWeight','bold',...
'String','START INVERSION','Tag','push1','CallBack',@main_inv)
%%%checkbox for instant RMS plot (if checked then time loss..)
uicontrol('Parent',panel1,'Style','checkbox','units','normalized',...
'Position',[0.05 0.54 0.4 0.03],'ForeGroundColor','r',...
'string','Instant RMS Plot','TooltipString',...
' this will cause time loss but supply visual control',...
'BackGroundColor',[.90 .90 .90],'FontWeight','bold',...
'Tag','checkbut1')

%%%%%
%%%%%LIST OF OUTPUT Graphics (plot by clicking)
panel3=uipanel('Parent',panel1,'Position',[0.05,.05,0.9,.45],...
'BackgroundColor',[.96 .96 .96],'Tag','pan3','BorderType','beveledin',...
'Title','Output Graphics & Plot Options','ForegroundColor','k','TitlePosition','centertop');
%%%text area for tic toc
uicontrol('Parent',panel3,'Style','text','units','normalized',...
'Position',[0.02 0.95 0.96 0.04],...
'FontWeight','normal','HorizontalAlignment','left','BackGroundColor',[.96 .96 .96],...
'ForeGroundColor','b','string','Elapsed time (sec.): 0','Tag','txttic')
%%% listbox of outs
uicontrol('Parent',panel3,'Style','listbox','units','normalized',...
'Position',[.02,.65,.96,.26],...
'String',{'Observed Anomaly';'Recalculated Anomaly';'Basement Surface';...
'Anomaly Difference';'RMS Variation'},'Tag','listb2',...
'FontWeight','bold','CallBack',@selectionoutput)
%%%%%%%Plotting Obtions
bg2=uibuttongroup(panel3,'units','normalized','Position',[0.01 0.35 0.979 0.25],...
'Title','Mapping Style','TitlePosition','lefttop','SelectionChangeFcn',@mapcontr);
%%%radiobuttons in button group
uicontrol(bg2,'Style','radiobutton','String','Filled Contour',...
'units','normalized','Position',[.05 .6 .45 .3],'Tag','rbut3');
uicontrol(bg2,'Style','radiobutton','String','Contour',...
'units','normalized', 'Position',[.55 .6 .40 .3],'Tag','rbut4');
uicontrol(bg2,'Style','radiobutton','String','3-D Filled Surface',...
'units','normalized', 'Position',[.05 .2 .5 .3],'Tag','rbut5');
uicontrol(bg2,'Style','radiobutton','String','3-D Surface',...
'units','normalized', 'Position',[.55 .2 .35 .3],'Tag','rbut6');
%%%%%%%%%%%%%%%%colormap choice
uicontrol('Parent',panel3,'Style','pushbutton','units','normalized',...
'Position',[0.4 0.27 0.58 0.05],...
'FontWeight','bold','BackGroundColor',[.90 .90 .90],...
'string','Open Colormap Editor','CallBack',@clrmaped)
%%%%%%%%%%%%%%%%%%%%%%
%%%% GET Anomaly PROFILE (+Cross Section of Basement Depth)  
uicontrol('Parent',panel3,'Style','pushbutton','units','normalized',...
'Position',[0.05 0.02 0.90 0.15],'ForeGroundColor','b',...
'FontWeight','bold',...
'String','GET A CROSS-SECTION','Tag','push2','CallBack',@crossec)

%%%%%uitable
uitable('Parent',fig1,'units','normalized','Position',[.2 .9 .76 .085],...
'RowName',{'RMS'},'Tag','tabl1','visible','off');
%%%%SET MENU ITEMS enable OFF at START
set(findall(panel1, '-property', 'enable'), 'enable', 'off')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Listbox for forward calc-mode
uicontrol('Parent',fig1,'Style','listbox','units','normalized',...
'Position',[.2,.94,.1,.05],...
'String',{'Depth-Model';'Gravity Anomaly'},'Tag','listbmodel',...
'visible','off','FontWeight','bold','CallBack',@selectionmodel)
%%%%%pushbutton save forward model
uicontrol('Parent',fig1,'Style','pushbutton','units','normalized',...
'Position',[.3,.94,.08,.03],...
'String','Save Grid','Tag','savlistbmodel',...
'visible','off','FontWeight','bold','CallBack',@savlistmodel)
%%%%%

%%%%%% Create an Axis 
axes('Parent',fig1,'units','normalized','Position',...
    [0.3 0.12 0.6 0.8],'Tag','axsig');
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%end of main figure 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAIN CORE CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main_inv(src,event)
%%%temporary plot and settings
plot(.5,.5,'k+','MarkerSize',10);
title('Evaluating...please wait');axis off;drawnow
set(findobj(gcf,'Tag','listb2'),'ForeGroundColor','k','value',5)
set(findobj(gcf,'Tag','push1'),'ForeGroundColor','r')
set(findobj(gcf,'Tag','txttic'),'ForeGroundColor','r')
set(findobj(gcf,'Tag','tabl1'),'Data',[],'visible','off');
drawnow
%%%%ask for maximum number of iteration
prompt = {'Maximum Nunber of Iteration'};
dlg_title = 'input';
num_lines = 1;
def = {'50'};
answer= inputdlg(prompt,dlg_title,num_lines,def);
itermax=str2num(cell2mat(answer));
%import active grid
recdat=load('gchin.mat');
nx=recdat.nx;
ny=recdat.ny;
g0=recdat.g0;
x=recdat.x;
y=recdat.y;
%square actual spacing interval in km...dx=dy; 
dx=str2num(get(findobj(gcf,'Tag','edactdx'),'string'));
dy=dx;
%density contrast
r0=str2num(get(findobj(gcf,'Tag','edro'),'string'));
%decay constant
lambda=str2num(get(findobj(gcf,'Tag','edlambda'),'string'));
%%ignoring any data node?
nig=ceil(str2num(get(findobj(gcf,'Tag','edign'),'string')));
%%%%
%order
n=6;
%iteration control
itcon=get(findobj(gcf,'Tag','rbut2'),'value');
% %%if itcon=1 iteration threshold, if 0 then Auto.
if itcon==1;
    criterio=str2num(get(findobj(gcf,'Tag','edcriterio'),'string')); 
    rms=ceil(criterio+1);
else
    criterio=0;
end
%%% Supply instant RMS plot?(if yes time loss for tic toc)
checkrms=get(findobj(gcf,'Tag','checkbut1'),'value');
%%%
const=41.9047*r0;
z=-1./lambda.*log(1-(lambda.*g0./const));
iter=0;
%%%%%inversion procedure starts (mode: iteration stop rms<threshold)
if itcon==1;
if r0<0    
tic
%rms=10;
while rms>criterio
iter=iter+1;
g=FW_Granser(z,r0,lambda,nx,ny,dx,dy,n);
dg=g0(nig+1:ny-nig,nig+1:nx-nig)-g(nig+1:ny-nig,nig+1:nx-nig);
nnxny=size(dg);
dg2=dg.^2;
rms=sqrt(sum(sum(dg2))./((nnxny(1)*nnxny(2))));
zold=z;%%store the actual z of the latest g
z=g0./g.*z;
r(iter)=rms;
if checkrms==1
plot(r(1:end),'-ks','MarkerFaceColor','b','MarkerSize',4)
drawnow
end
if iter==itermax;break;end
end
time=toc;
z=zold;
else
    msgbox('density contrast ?')
end
end
%%%%%Convergence Threshold..iteration procedure ends

%%%%%inversion procedure starts (mode: divergence..iteration stops when Rms(i+1)>rms(i))
if itcon==0;
if r0<0
r=1000*ones(1,itermax);
%%%%%%%%%%%%%%input dialog box for maximum iteration
tic
for j=2:itermax
g=FW_Granser(z,r0,lambda,nx,ny,dx,dy,n);
dg=g0(nig+1:ny-nig,nig+1:nx-nig)-g(nig+1:ny-nig,nig+1:nx-nig);
nnxny=size(dg);
dg2=dg.^2;
rms=sqrt(sum(sum(dg2))./((nnxny(1)*nnxny(2))));
r(j)=rms;
if r(j)>r(j-1) || j==itermax
    if j>2;g=gold;z=zold;dg=dgold;end
    break
else
gold=g;
zold=z;
dgold=dg;
z=g0./g.*z;
if checkrms==1
plot(r(2:j),'-ks','MarkerFaceColor','b','MarkerSize',4)
drawnow
end
end
end
time=toc;
findex=find(r==1000);r(findex)=[];
else
    msgbox('density contrast ?')
end
end
%%%%% divergence mode iteration procedure ends
%%%save outputs to a temporary mat file 
g=g(nig+1:ny-nig,nig+1:nx-nig);
g0=g0(nig+1:ny-nig,nig+1:nx-nig);
z=z(nig+1:ny-nig,nig+1:nx-nig);
x=x(nig+1:ny-nig,nig+1:nx-nig);
y=y(nig+1:ny-nig,nig+1:nx-nig);
save('gchout.mat','g','g0','dg','z','r','x','y','dx','dy','n','lambda','r0','criterio');
%settings;
set(findobj(gcf,'Tag','txttic'),'String',['Elapsed time (sec.):' ...
num2str(time,'%10.2f')])
set(findobj(gcf,'Tag','listb2'),'ForeGroundColor','b','value',5)
set(findobj(gcf,'Tag','push1'),'ForeGroundColor','k')
set(findobj(gcf,'Tag','txttic'),'ForeGroundColor','k')
%%%run listbox for plotting
selectionoutput
set(findobj(gcf,'Tag','savbut'),'enable','on');
set(findobj(gcf,'Tag','expbut'),'enable','on');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Forward Granser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g=FW_Granser(z,r0,lambda,nx,ny,dx,dy,n)
z0=(max(max(z))-min(min(z)))/2;
%%  extend the model area
nx0=nx; ny0=ny;
z(1,nx+floor(nx/2))=0;
z(ny+floor(ny/2),1)=0;
z=rot90(rot90(z));
z(1,nx+2*floor(nx/2))=0;
z(ny+2*floor(ny/2),1)=0;
z=rot90(rot90(z));
if (mod(nx,2)~=0); nx=nx-1; z(:,end)=[];end
if (mod(ny,2)~=0); ny=ny-1; z(end,:)=[];end

%%
nxm=2*nx; nym=2*ny;
h=z-z0;
dkx= 2.*pi./((nxm-1).*dx);
dky= 2.*pi./((nym-1).*dy);
nyqx= (nxm/2)+1;
nyqy= (nym/2)+1; 
for j=1:nxm
       if j <= nyqx
         kx(j)=(j-1)*dkx;
       else
         kx(j)=(j-(nxm+1))*dkx;
       end 
for i=1:nym
       if i <= nyqy
         ky(i)=(i-1)*dky;
       else
         ky(i)=(i-(nym+1))*dky;
       end 
       k(i,j)=sqrt(kx(j).^2+ky(i).^2);
end
end
hs1=2*pi*20/3.*r0.*exp(-lambda.*z0);
hs2=exp(-abs(k).*z0)./(lambda+abs(k));
tongF=fft2(1-exp(-lambda.*h));
for m=1:n;
     tongF=tongF-((-k).^(m))./(factorial(m)).*fft2(exp(-lambda.*h).*h.^m); %
end;
Fg=hs2.*tongF;
g0=(ifft2(Fg));
r_g0=real(g0);
g1=r_g0.*hs1;
g1=g1(ny/2+1:ny/2+ny0,nx/2+1:nx/2+nx0);
g=g1+2.*pi.*20/3.*r0./lambda.*(1-exp(-lambda.*z0));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%ADDITIONAL FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%CROSS SECTION
function crossec(src,event)
set(findobj(gcf,'Tag','rbut3'),'value',1)
set(findobj(gcf,'Tag','listb2'),'value',3);
selectionoutput
drawnow
recin=load('gchin.mat');
recout=load('gchout.mat');
[xxx,yyy]=getline;
if length(xxx)>1
xp(1)=xxx(1);yp(1)=yyy(1);
xp(2)=xxx(2);yp(2)=yyy(2);
point1 = [xp(1);yp(1)];ppx=[xp(1) xp(2)];
point2 = [xp(2);yp(2)];ppy=[yp(1) yp(2)];
proflen=sqrt((abs(xp(1)-xp(2)))^2+(abs(yp(1)-yp(2)))^2);
%%%%%divide to point number
prompt = {['Divide section to N point (section length is ' num2str(proflen)]};
dlg_title = 'setting';
num_lines = 1;
def = {'100'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
npoint=str2num((cell2mat(answer(1))));
%%%%%%%%%%%%%%%%%%%%%%%%%%
np=npoint -1;
ndx=1/np;
t = 0:ndx:1;
xi = point1(1)*t+point2(1)*(1-t);
yi = point1(2)*t+point2(2)*(1-t);
zi=interp2(recin.x,recin.y,recin.g0,xi,yi);
zi=fliplr(zi);
zii=interp2(recout.x,recout.y,recout.g,xi,yi);
zii=fliplr(zii);
ziii=interp2(recout.x,recout.y,-recout.z,xi,yi);
ziii=fliplr(ziii);
hold on
nn=length(xi)-1;dx=proflen/nn;
line(xi,yi,'linewidth',4,'color','k');
hold off
figure('MenuBar','none','Name',...
'GCH_gravity','NumberTitle','off','ReSize','off',...
'Color','w',...
'units','normalized','outerposition',[0 0.5 .45 .45]);
xpro=0:dx:proflen;
plot(xpro,zi,'ro','MarkerFaceColor','r','MarkerSize',2);hold on
plot(xpro,zii,'-k','LineWidth',2);hold off;grid on
xlabel('Distance');ylabel('mGal');legend('Observed','Calculated');
title('Cross-Sections of Observed and Recalculated Anomalies')
xlim([min(xpro) max(xpro)])
figure('MenuBar','none','Name',...
'GCH_depth','NumberTitle','off','ReSize','off',...
'Color','w',...
'units','normalized','outerposition',[0 0.03 .45 .45]);
plot(xpro,ziii,'-k','LineWidth',2);grid on
karex1=[xpro xpro(end) xpro(1) xpro(1)];
karey1=[ziii 0 0 ziii(1)];
karex2=[xpro xpro(end) xpro(1) xpro(1)];
karey2=[ziii min(ziii) min(ziii) ziii(1)];
patch(karex1,karey1,[.8 .8 .8],'FaceAlpha',.5);
patch(karex2,karey2,[.4 .4 .4],'FaceAlpha',.5);
xlabel('Distance');ylabel('Depth');
xlim([min(xpro) max(xpro)])
ylim([min(ziii) 0])
pbaspect([1 .3 1])
title('Cross-Section of Basement Depth')

end
end

%%%%saving results to a mat file
function sav2data(src,event)
k=get(findobj(gcf,'Tag','mainfig'),'name');
[PATHSTR,NAME,EXT] = fileparts(k);
[filename, pathname] = uiputfile([NAME '_gch.mat'], 'Save Output data');
kk=[pathname filename];
if ischar(kk)
copyfile('gchout.mat',kk)
msgbox('Data Stored as mat file...')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%saving results to a mat file
function sav1data(src,event)
k=get(findobj(gcf,'Tag','mainfig'),'name');
[PATHSTR,NAME,EXT] = fileparts(k);
[filename, pathname] = uiputfile([NAME '_rec.grd'], 'Save Output data to *.grd');
kk=[pathname filename];
if ischar(kk)
rec=load('gchout.mat');
xmin=min(min(rec.x));ymin=min(min(rec.y));
xmax=max(max(rec.x));ymax=max(max(rec.y));
grdout(rec.g,xmin,xmax,ymin,ymax,[filename(1:end-4) '_gch_calc-g.grd'])
grdout(-rec.z,xmin,xmax,ymin,ymax,[filename(1:end-4) '_gch_calc-z.grd'])
grdout(rec.dg,xmin,xmax,ymin,ymax,[filename(1:end-4) '_gch_calc-diffg.grd'])
rms=rec.r;
rmsn=[(1:numel(rms))' rms'];
save([filename(1:end-4) '_gch-rms.dat'],'rmsn','-ascii')
msgbox('DATA FILES ARE STORED with extension *calc-g, calc-z and *rms')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%radiobuttons behavior
function radbehav(src,event)
val2=get(findobj(gcf,'Tag','rbut2'),'value');
if val2==1;set(findobj(gcf,'Tag','edcriterio'),'enable','on');end
if val2==0;set(findobj(gcf,'Tag','edcriterio'),'enable','off');end
end

function mapcontr(src,event)
selectionoutput
end

%%%%Get output selection for plotting
function selectionoutput(src,event)
rf2=exist('gchout.mat','file');
if rf2==2 
recout=load('gchout.mat');
val=get(findobj(gcf,'Tag','listb2'),'value');
%%%get plot obtions
v1=get(findobj(gcf,'Tag','rbut3'),'value');%%contourf
v2=get(findobj(gcf,'Tag','rbut4'),'value');%%contour
v3=get(findobj(gcf,'Tag','rbut5'),'value');%%colored surf
v4=get(findobj(gcf,'Tag','rbut6'),'value');%%nocolor surf
stylemap=[v1 v2 v3 v4];index=find(stylemap==1);

switch val
    case 1
        mapper(recout.x,recout.y,recout.g0,'mGal','Observed Anomaly',index)
        set(findobj(gcf,'Tag','tabl1'),'Data',[],'visible','off');
        case 2
        mapper(recout.x,recout.y,recout.g,'mGal','Inverted Anomaly',index)
        set(findobj(gcf,'Tag','tabl1'),'Data',[],'visible','off');
        case 3
        mapper(recout.x,recout.y,-recout.z,'km','Basement Surface',index)
        set(findobj(gcf,'Tag','tabl1'),'Data',[],'visible','off');
        case 4
        mapper(recout.x,recout.y,recout.dg,'mGal','Differentation',index)
        set(findobj(gcf,'Tag','tabl1'),'Data',[],'visible','off');
        case 5
        plot(recout.r,'-ks','MarkerFaceColor','r','MarkerSize',5)
        grid on;xlabel('Iteration Number');ylabel('RMS');
        rindex=find(recout.r==min(recout.r));
        title(['RMS Variation: optimum rms=' num2str(min(recout.r)) ...
            '  at iteration step:' num2str(rindex)])
        set(findobj(gcf,'Tag','tabl1'),'Data',recout.r,'visible','on');
end
else
    msgbox('Starting of inversion required...')
end
end

function mapper(x,y,matrix,unitt,tit,index)

if index==1;
    contourf(x,y,matrix,30);shading flat;rotate3d off;end
if index==2;
    contour(x,y,matrix,30);rotate3d off;end
if index==3;
surf(x,y,matrix,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud')
camlight right
rotate3d on
end

if index==4;
mesh(x,y,matrix)
rotate3d on
end
h=colorbar('eastoutside');title(h,unitt);
xlabel('Distance X');ylabel('Distance Y');title(tit)
grid on
rf1=load('gchin.mat');
xlim([min(min(rf1.x)) max(max(rf1.x))])
ylim([min(min(rf1.y)) max(max(rf1.y))])
axis equal
axis tight

%%%%%messagebox
if index==3 || index==4
hm=uicontrol('Parent',gcf,'Style','text','units','normalized',...
'Position',[.5 .5 .1 .05],'string','Rotate On','BackGroundColor','k',...
'ForeGroundColor','g','FontWeight','bold','FontSize',16);
pause(2);delete(hm);
end
end

function clrmaped(src,event)
colormapeditor
end

%%%%%action menu 1 IMPORT a *.grd file (created in format surfer 6 text grid)
function impgrd(src,event);
clc;
[filename, pathname] = uigetfile('*.grd', 'Import Golden Software grid (*.grd)');
k=[pathname filename];
if ischar(k)
fidc=fopen(k);
header= fread(fidc,4,'*char' )';
fclose(fidc);
c1=strcmp(header,'DSAA');
c2=strcmp(header,'DSRB');
sumc=sum([c1 c2]);
if sumc>0
%%delete any temporary in and out files of GCHinv..
if exist('gchin.mat','file')==2;delete('gchin.mat');end
if exist('gchout.mat','file')==2;delete('gchout.mat');end
set(findobj(gcf,'Tag','mainfig'),'Name',k)
switch c1
    case 1
[g0,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(k);
    case 0
[g0,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd7bin(k);        
end
list1w(g0,nx,ny,xmin,xmax,ymin,ymax,dx,dy,k);
save('gchin.mat','g0','x','y','nx','ny','xmin','xmax','ymin','ymax','dx','dy');
set(findobj(gcf,'Tag','listb2'),'ForeGroundColor','k')
set(findobj(gcf,'Tag','push1'),'ForeGroundColor','b')
set(findobj(gcf,'Tag','txttic'),'ForeGroundColor','b')
set(findobj(gcf,'Tag','tabl1'),'Data',[],'visible','off');
plot(.5,.5,'k+')
axis off
drawnow
mapper(x,y,g0,'mGal','Observed Anomaly',1)
set(findall(findobj(gcf,'Tag','pan1'), '-property', 'enable'), 'enable', 'on')
set(findobj(gcf,'Tag','listbmodel'),'value',1,'visible','off');
set(findobj(gcf,'Tag','savlistbmodel'),'visible','off');
set(findobj(gcf,'Tag','savbut'),'enable','off');
set(findobj(gcf,'Tag','expbut'),'enable','off');
else
msgbox('File format not supported..Load Surfer 6 text grid / Surfer 7 Binary Grid')
end
end
end
%%%%%action menu 1 IMPORT a XYZ (*.dat) equal spaced grid file 
function impgrdat(src,event);
clc;
[filename, pathname] = uigetfile('*.dat', 'Import XYZ grid data(*.dat)');
k=[pathname filename];
if ischar(k)
if exist('gchin.mat','file')==2;delete('gchin.mat');end
if exist('gchout.mat','file')==2;delete('gchout.mat');end    
set(findobj(gcf,'Tag','mainfig'),'Name',k)
[g0,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrdatfile_Callback(k);
list1w(g0,nx,ny,xmin,xmax,ymin,ymax,dx,dy,k)
save('gchin.mat','g0','x','y','nx','ny','xmin','xmax','ymin','ymax','dx','dy');
set(findobj(gcf,'Tag','listb2'),'ForeGroundColor','k')
set(findobj(gcf,'Tag','push1'),'ForeGroundColor','b')
set(findobj(gcf,'Tag','txttic'),'ForeGroundColor','b')
set(findobj(gcf,'Tag','tabl1'),'Data',[],'visible','off');
plot(.5,.5,'k+')
axis off
drawnow
mapper(x,y,g0,'mGal','Observed Anomaly',1)
set(findall(findobj(gcf,'Tag','pan1'), '-property', 'enable'), 'enable', 'on')
set(findobj(gcf,'Tag','listbmodel'),'value',1,'visible','off');
set(findobj(gcf,'Tag','savlistbmodel'),'visible','off');
set(findobj(gcf,'Tag','savbut'),'enable','off');
set(findobj(gcf,'Tag','expbut'),'enable','off');
end
end

%%%%%read Surfer 6 text grid(*.grd)
function [g0,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(k)
surfergrd=fopen(k,'r'); % Open *.grid file
dsaa=fgetl(surfergrd);  % Header
% Get the map dimension [NX: East NY: North];
datasize=str2num(fgetl(surfergrd)); nx=datasize(1); ny=datasize(2);
% Map limits: xmin, xmax, ymin ymax
 xcoor=str2num(fgetl(surfergrd)); xmin=xcoor(1); xmax=xcoor(2);
 ycoor=str2num(fgetl(surfergrd)); ymin=ycoor(1); ymax=ycoor(2);
% check intervals in x and y direction 
dx=(xmax-xmin)/(nx-1);
dy=(ymax-ymin)/(ny-1);
% data limits
anom=str2num(fgetl(surfergrd)); g0min=anom(1); g0max=anom(2);
% data matrix 
[g0,numb] = fscanf(surfergrd, '%f', [nx,ny]);
g0=g0'; % Traspose matrix
fclose(surfergrd);
% map coordinate matrix
[x,y]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read Surfer 7 Binary grid
function [g0,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy] = lodgrd7bin(filename)
fid= fopen(filename);
fread(fid,4,'*char' )';
fread(fid,1,'uint32');fread(fid,1,'uint32');
fread(fid,4,'*char' )';fread(fid,1,'uint32');
ny= fread(fid,1,'uint32'); nx= fread(fid,1,'uint32');
xmin= fread(fid,1,'double'); ymin= fread(fid,1,'double');
dx= fread(fid,1,'double'); dy= fread(fid,1,'double');
fread(fid,1,'double');fread(fid,1,'double');
fread(fid,1,'double');
parm= fread(fid,1,'double');
fread(fid,4,'*char' )';
nn= fread(fid,1,'uint32');
if ny*nx ~= nn/8 ; error('error') ;end
g0= nan(nx,ny);
g0(1:end) = fread(fid,numel(g0),'double');
g0=g0';
fclose(fid);
g0(g0==parm) = nan;
xv = xmin + (0:nx-1)*dx;
yv = ymin + (0:ny-1)*dy;
[x,y]=meshgrid(xv,yv);
xmax=xv(end);
ymax=yv(end);
end
%%%%%%%%%%read XYZ gridded data (*.dat)
function [g0,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrdatfile_Callback(k)
data=load(k);
x=data(:,1);y=data(:,2);g0=data(:,3); 
unix=unique(x);xmin=unix(1);xmax=unix(end);dx=unix(2)-unix(1);
uniy=unique(y);ymin=uniy(1);ymax=uniy(end);dy=uniy(2)-uniy(1);
nx=numel(unix);
ny=numel(uniy);
if x(2)==x(1);
g0=reshape(g0,ny,nx);
x=reshape(x,ny,nx);
y=reshape(y,ny,nx);
else
g0=reshape(g0,nx,ny);
x=reshape(x,nx,ny);
y=reshape(y,nx,ny);
g0=g0';x=x';y=y';
end
end

%%%%%%%%%%% Function for output of a GRID
function grdout(matrix,xmin,xmax,ymin,ymax,namefile)
%Get grid dimensions
aux=size(matrix);
nx=aux(2);ny=aux(1);
grdfile=fopen(namefile,'w');                % Open file
fprintf(grdfile,'%c','DSAA');               % Header code
fprintf(grdfile,'\n %i %i',[nx ny]);        % Grid size
fprintf(grdfile,'\n %f %f',[xmin xmax]);    % X limits
fprintf(grdfile,'\n %f %f',[ymin ymax]);    % Y limits
fprintf(grdfile,'\n %f %f',[min(min(matrix)) max(max(matrix))]); % Z limits
fprintf(grdfile,'\n');
for jj=1:ny                                 % Write matrix
    for ii=1:nx
        fprintf(grdfile,'%g %c',matrix(jj,ii),' ');
    end
    fprintf(grdfile,'\n');
end
fclose(grdfile);
end

%%%%%%%%%%%%%%%%%Get Grid Info and write it to listbox1
function list1w(g0,nx,ny,xmin,xmax,ymin,ymax,dx,dy,k)
set(findobj(gcf,'Tag','listb1'),'value',1);
s1=['Source : ' k];
s2=['NX : ' num2str(nx)];
s3=['NY : ' num2str(ny)];
s4=['Grid Spacing dx (in map units): ' num2str(dx)];
s5=['Grid Spacing dy (in map units): ' num2str(dy)];
s6=['Xmin - Xmax:   ' num2str(xmin) '  /   ' num2str(xmax)];
s7=['Ymin - Ymax:   ' num2str(ymin) '  /   ' num2str(ymax)];
s8=['Zmin - Zmax:   ' num2str(min(min(g0))) '  /  ' num2str(max(max(g0)))];
str={s1;s2;s3;s4;s5;s6;s7;s8};
set(findobj(gcf,'Tag','listb1'),'string',str);
set(findobj(gcf,'Tag','edactdx'),'string',num2str(dx));
end

%%%%%%Export Figure as..
function expfig(src,event)
val=get(findobj(gcf,'Tag','listb2'),'value');
hax = findall(gcf,'type','axes'); % Find the axes object in the GUI
cm = get(gcf,'Colormap');
f1 = figure('Color','w','Name','Set & Save Figure','NumberTitle','off') ;
copyobj(hax(end),f1); % Copy axes object h into figure f1
if val<5;set(gcf,'colormap',cm);colorbar;end
end
%%%%
%%%%%Forward calculation of gravity from a defined depth grid;
function forwardc(src,event)
clc;
[filename, pathname] = uigetfile('*.grd', 'Import Depth Grid (*.grd)');
k=[pathname filename];
if ischar(k)
fidc=fopen(k);
header= fread(fidc,4,'*char' )';
fclose(fidc);
c1=strcmp(header,'DSAA');
c2=strcmp(header,'DSRB');
sumc=sum([c1 c2]);
if sumc>0
switch c1
    case 1
[z,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(k);
    case 0
[z,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd7bin(k);        
end
prompt = {'Density constrast (negative!)','Decay constant'};
dlg_title = 'Model Parameters';
num_lines = 1;
def = {'-0.57','0.25'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r0=str2num(cell2mat(answer(1)));
lambda=str2num(cell2mat(answer(2)));
n=6;
z=abs(z); %Depth is positive downward
g=FW_Granser(z,r0,lambda,nx,ny,dx,dy,n);
save('gchz2g.mat','z','g','r0','lambda','x','y','nx','ny','xmin','xmax',...
    'ymin','ymax','dx','dy','filename');
plot(.5,.5,'k+')
axis off
drawnow
mappermod(x,y,-abs(z),'km','Depth-Model')
xlim([xmin xmax])
ylim([ymin ymax])
set(findall(findobj(gcf,'Tag','pan1'), '-property', 'enable'), 'enable', 'off')
set(findobj(gcf,'Tag','tabl1'),'visible','off');
set(findobj(gcf,'Tag','listbmodel'),'value',1,'visible','on','enable','on');
set(findobj(gcf,'Tag','savlistbmodel'),'visible','on','enable','on');
set(findobj(gcf,'Tag','savbut'),'enable','off');
set(findobj(gcf,'Tag','expbut'),'enable','off');
else
msgbox('File format not supported..Load Surfer 6 text grid / Surfer 7 Binary Grid')
end
end
end

function mappermod (x,y,matrix,unitt,tit)
contourf(x,y,matrix,30);shading flat;rotate3d off
h=colorbar('eastoutside');title(h,unitt);
xlabel('Distance X');ylabel('Distance Y');title(tit)
grid on
axis equal
axis tight
end

function selectionmodel(scr,event)
recor=load('gchz2g.mat');
val=get(findobj(gcf,'Tag','listbmodel'),'value');
switch val
    case 1
        mappermod(recor.x,recor.y,-(abs(recor.z)),'km','Depth-Model')
        set(findobj(gcf,'Tag','savlistbmodel'),'visible','on','enable','off');
    case 2
    tit=['Gravity-Model : density contrast=' num2str(recor.r0) ...
        '  decay constant=' num2str(recor.lambda)];
        mappermod(recor.x,recor.y,recor.g,'mGal',tit)
        set(findobj(gcf,'Tag','savlistbmodel'),'visible','on','enable','on');
end
xlim([recor.xmin recor.xmax])
ylim([recor.ymin recor.ymax])
end

function savlistmodel(src,event)
rec=load('gchz2g.mat');
[filename, pathname] = uiputfile([rec.filename(1:end-4) '_modg.grd'], ...
    'Save Output data to *.grd');
kk=[pathname filename];
if ischar(kk)
xmin=min(min(rec.x));ymin=min(min(rec.y));
xmax=max(max(rec.x));ymax=max(max(rec.y));
grdout(rec.g,xmin,xmax,ymin,ymax,kk);
set(findobj(gcf,'Tag','listbmodel'),'value',1,'enable','off');
set(findobj(gcf,'Tag','savlistbmodel'),'enable','off');
msgbox('Gravity model has been stored...')
end
end
