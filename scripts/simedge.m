function simedge(year,meta)
res=1000;
days=365+mod(year,4);
mesh=zeros(days, floor(res/2));
sunarray=generate_sun_array(year,res,meta);
for i=1:days
    for j=floor(res/4):floor(3*res/4)
        offset=floor(res/4)-1;
        if (sunarray(i,j,2)>15&&sunarray(i,j,2)<30)
            mesh(i,j-offset)=1;
        end
    end
end
edge=diff(mesh);
[row,col]=find(edge);
col=col/1000*48;
plot(row,col)







            

%object description
%object az>lhs az<rhs zenith any square object
%eliptical object 
%simulates an object edge and predicts the objects shadowing edge 
