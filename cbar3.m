function cbar3(A)

h=bar3(A); % Make the standard plot....
Max_A=max(A(:)); % Tweak with the colormap/colorbar...
Min_A=min(A(:));
caxis([Min_A Max_A]);
colormap(jet);
colorbar;


for j=1:size(h,2) % For each handle returned from bar3, NOTE size(h,2)=size(A,2);
 
  CData=get(h(j),'CData'); % Get the color data
  ZData=get(h(j),'ZData'); % Get the height data
 
  Size_A=size(A,1);
 
  for k=1:Size_A % Loop over each "set" of bars....
    magic_ind=6*(k-1)+1;
    Z_val=ZData(magic_ind+1,2); % Determine the Z-value for each bar
    CData(magic_ind:k*6,:)=Z_val; % Set the color of the bar based on the Z-value
  end
 
  set(h(j),'CData',CData); % Change the CData handle to reflect the different colors
 
end


return

