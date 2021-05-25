function resampled_polyg=resample_polyg(polyg, npts)
% Juan Camilo Gutierrez U
% 201710009014
% 11/04/2021

% pts = (2xN) or (3xN) polygon to be resampled
% npts = number of point per polygon segment

sides = size(polyg,2); 

resampled_polyg = [];
for i= 1:sides
    a = polyg(:,i);
    if i == sides
        b = polyg(:,1);
    else
        b = polyg(:,i+1);   
    end    
    resampled_line = [];
    for j = 0:npts-1
        delta = ((b-a)/npts)*j;
        npoint = polyg(:,i)+delta;
        resampled_line = [resampled_line npoint];
    end
    resampled_polyg = [resampled_polyg resampled_line];
end


end