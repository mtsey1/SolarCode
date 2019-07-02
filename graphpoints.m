function [xx,yy,indx,sect]=graphpoints(varargin)
% graphpoints  Get coordinates and indices of points from a graph
%  2016-10-12  Matlab2016  
%  Copyright (C) 2016, W.Whiten (personal W.Whiten@uq.edu.au) BSD license
%
% [xx,yy,indx,sect]=graphpoints()            % uses current figure
% [xx,yy,indx]=graphpoints(x,y,options)      % generates figure
% [xx,yy,indx]=graphpoints(x,y,ind,options)  % selected points on figure
%  x     X coordinates of points to plot
%  y     Y coordinates    Give x,y or x,y,ind or omit all three
%  ind   Indices of points to plot  (optional)
%  options Name value pairs, or struct, see optndfts
%          .x       X coordinates  (default from current figure)
%          .y       Y coordinates 
%          .ind     Indices of points (default 1:n)
%          .symbol  Code for point to plot (default '*r')
%          .prompts Option to include screen prompts (default true)
%
%  xx    X coordinates of selected values
%  yy    Y coordinates of selected values
%  indx  Indices of points selected, within section if miultiple sections
%  sect  Number of section when graph contains multiple sections
%
% Plots x,y or or x(ind),y(ind) or uses points from current figure.
% Select points from graph with the cross lines, press return to 
%  terminate point collection.
% When prompted "Exit graphpoints (y):" graph can be enlarged or rescaled, 
%  and any answer except initial y will return to collecting more points 
%  from graph. Set option 'prompts' to false to skip keyboard prompts.
% If graph has multiple sections (h=gcf;length(h.CurrentAxes.Children)>1) 
%  indx is count in within section, and sect is section number. Note
%  section number one is last section added to graph.

optn=optndfts(varargin,{'x','y','ind'},'x',[],'y',[],'ind',[],  ...
    'symbol','+','text',true);
x=optn.x;
y=optn.y;
ind=optn.ind;
if(islogical(ind))
    ind=find(ind);
end

if(optn.text)
    disp(' ')
    disp('Place curser near each point, and click, ')
    disp('use "return" to exit graph')
    pause(2);
end

if(isempty(x))
    
    % get point coordinates from current figure
    h=gcf;
    h=h.CurrentAxes.Children;
    nh=length(h);
    
    nx=0;
    for i=1:nh
        nx=nx+length(h(i).XData);
    end
    x=zeros(nx,1);
    y=zeros(nx,1);
    sectcnt=zeros(nh,1);
        
    cnt1=1;
    for i=1:nh
        cnt2=cnt1+length(h(i).XData);
        x(cnt1:cnt2-1)=h(i).XData;
        y(cnt1:cnt2-1)=h(i).YData;
        sectcnt(i)=cnt1;
        cnt1=cnt2;
    end
    if(isempty(ind))
        ind=1:length(x);
    end
else
    
    % create new figure
    if(isempty(ind))
        ind=1:length(x);
    end
    hold off
    plot(x(ind),y(ind),optn.symbol);
    sectcnt=1;
end

% cycles of get from graph
j=0;
mj=100;
xx=zeros(mj,1);
yy=zeros(mj,1);
indx=zeros(mj,1);
sect=zeros(mj,1);
while(true)
    
    % get points from graph
    pnts=ginput;
    
    % peocess point selected
    for i=1:size(pnts,1)
        j=j+1;
        if(j>mj)    % expand result vectors
            mj=round(mj*1.5);
            xx(mj)=0;
            yy(mj)=0;
            indx(my)=0;
            sect(mj)=0;
        end
        
        % locate nearest point and save
        d=(x(ind(:))-pnts(i,1)).^2+(y(ind(:))-pnts(i,2)).^2;
        [~,i1]=min(d);
        xx(j)=x(ind(i1));
        yy(j)=y(ind(i1));
        indx(j)=ind(i1);
        
        sect(j)=sum(sectcnt<=indx(j));
        indx(j)=indx(j)-sectcnt(sect(j))+1;
    end
    
    % test if more input requested
    if(optn.text)
        str=input('Exit graphpoints (y):','s');
        if(strncmpi(str,'y',1))
            break
        end
    else
        break
    end
end

xx=xx(1:j,1);
yy=yy(1:j,1);
indx=indx(1:j,1);
sect=sect(1:j,1);

return
end
