function gen_start_fit(data_no_vamp,sunPos)

%%reshape to make sure there's no singleton dimensions
data_no_vamp=squeeze(data_no_vamp);
%%create an array that will contain our start/stop contours
start_stop = zeros(size(data_no_vamp,1),2);

for i = 1:size(data_no_vamp,1)
    check_start = 0;
    check_stop = 0;
    for j = 2:size(data_no_vamp,2)
        if data_no_vamp(i,j-1)<0&&data_no_vamp(i,j)>0&&check_start==0
           check_start = 1;
           start_stop(i,1)=j;
        end
            
        if data_no_vamp(i,j-1)>0&&data_no_vamp(i,j)<0&&check_stop==0
           check_stop = 1;
           start_stop(i,2)=j;
        end
    end
end


%%separate start/stop out into individual arrays, we'll need this to remove
%%sunrise and sunset times from each vector
start_cont = start_stop(:,1);
end_cont = start_stop(:,2);

%%need to remove zeroes, and remove times that coincide or preempt sunrise
    %%remove zeroes:
    %%we remove the zeroes, because they dont record any generation
start_cont = start_cont(start_cont>0);
end_cont = end_cont(end_cont>0);

    %%remove sunrise/sunset values




%%remove the times where start or stop coincides wit



end