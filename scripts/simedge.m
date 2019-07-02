function simedge(year,meta)
res=1000;
days=365+mod(year,4);
mesh=zeros(days, res);
sunarray=generate_sun_array(year,res,meta);
%might need to shift 0deg to south
sunnarray(:,:,2)=mod(sunarray(:,:,2)'+180,360);
leftheight=80;
rightheight=85;
leftedge=170;
rightedge=190;
centerheight=70;
for i=1:days
    for j=1:res
        if sunarray(i,j,3)>90 
            mesh=1;
        elseif (sunarray(i,j,3)< min(leftheight,rightheight))
            mesh(i,j)=1;
        elseif leftheight>rightheight
                if sunarray(i,j,2)>leftedge && sunnarray(i,j,3)<leftheight
                    mesh(i,j)=1;
                end
        elseif sunarray(i,j,2)<rightedge && sunnarray(i,j,3)<rightheight
                    mesh(i,j)=1;
        elseif sunarray(1,j,2)>rightedge && sunarray(1,j,2)<leftedge && sunnarray(i,j,3)>centerheight 
            mesh(i,j)=1;
        end 
    end
end
%itterative refining of paramaters


edge=diff(mesh);
[row,col]=find(edge);
col=col/1000*48;
plot(row,col)
end
