clear all;
clc;

load('pointTargetData.mat');

data = veraStrct.data(80:end,:,:);
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = 0.5*(speed/fs); 

for ii = 1:max(size(data))
    time_array_all(ii) = ii/fs;
end

for cc = 1:128
for bb = 1:128
    time_array(:,bb,cc) = time_array_all;
end
end

channel = [[-63.5:1:63.5]];

%2 parallel beams
data = veraStrct.data(80:end,:,1:2:128);
[rows_d col_d z_d] = size(data);

for beam = 1:z_d
    
for jj = 1:max(size(data)) %jj=row
    
depth = jj*pixel_size_through_depth; %m

data_matrix = data;
[rows_data_matrix col_data_matrix z_data_matrix] = size(data_matrix);

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    lat(ii) = 0.1953e-3*channel(ii);
    d_1(ii) = ((xe(ii)-((-1/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_1(ii) = d_1(ii)/speed;
    d_2(ii) = ((xe(ii)-((1/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_2(ii) = d_2(ii)/speed;
end

delay_matrix_1(jj,:,beam) = time_to_point_1; %delays
delay_matrix_2(jj,:,beam) = time_to_point_2;

end

for aa = 1:128
    delayed_channel_1(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_1(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_2(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_2(1:rows_data_matrix,aa,beam),'linear');
end


end
ax = [1:rows_data_matrix]*pixel_size_through_depth;

[delay_rows delay_col delay_beams] = size(delayed_channel_1);
delayed_channel = zeros(rows_data_matrix,128,128);
odd = 1;
even = 2;
for ind = 1:delay_beams %odd...delayed_channel_1
    delayed_channel(:,:,odd) = delayed_channel_1(:,:,ind);
    odd = odd+2;
    delayed_channel(:,:,even) = delayed_channel_2(:,:,ind);
    even = even+2;
end

for ll = 1:numel(delayed_channel)
    if isnan(delayed_channel(ll))==1
        delayed_channel(ll) = 0;
    end
end

summed_channels = sum(delayed_channel,2);
figure;
imagesc(lat, ax,20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
axis image;
title('2 parallel beams, (pointTargetData.mat)');


% %4 beams
% data = veraStrct.data(80:end,:,1:4:128);
% [rows_d col_d z_d] = size(data);

% for beam = 1:z_d
%     
% for jj = 1:max(size(data)) %jj=row
%     
% depth = jj*pixel_size_through_depth; %m
% 
% data_matrix = data;
% [rows_data_matrix col_data_matrix z_data_matrix] = size(data_matrix);
% 
% pitch = 0.1953e-3;
% xf = pitch*[-3/4 -1/4 1/4 3/4];
% for ll = 1:length(xf)
% for ii = 1:(length(channel))
%     xe(ii) = 0.1953e-3*abs(channel(ii)); 
%     d(ii) = ((xe(ii)-xf(ll))^2+depth^2)^0.5 + depth;
%     time_to_point(ii) = d(ii)/speed;
% end
% 
% delay_matrix(jj,:,beam) = time_to_point; %delays
% 
% cell_delays{ll,1} = delay_matrix;
% 
% end
% 
% end
% 
% end
% 
% n = 1;
% for beam_ind = 1:beam
%     delay_matrix = cell_delays{1,1};
%     for aa = 1:128
%         delayed_channel(1:rows_data_matrix,aa,n) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix(1:rows_data_matrix,aa,beam),'linear');
%     end
%     n=n+1;
%     delay_matrix = cell_delays{2,1};
%     for aa = 1:128
%         delayed_channel(1:rows_data_matrix,aa,n) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix(1:rows_data_matrix,aa,beam),'linear');
%     end
%     n=n+1;
%     delay_matrix = cell_delays{3,1};
%     for aa = 1:128
%         delayed_channel(1:rows_data_matrix,aa,n) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix(1:rows_data_matrix,aa,beam),'linear');
%     end
%     n=n+1;
%     delay_matrix = cell_delays{4,1};
%     for aa = 1:128
%         delayed_channel(1:rows_data_matrix,aa,n) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix(1:rows_data_matrix,aa,beam),'linear');
%     end
%     n=n+1;
% end
% 
% summed_channels = sum(delayed_channel,2);
% figure;
% imagesc(20*log10(abs(hilbert(summed_channels(:,:)))));
% colormap('gray');
% title('4 parallel beams, (pointTargetData.mat)');

% % load('anecoicCystData.mat');
% % 
% % data = veraStrct.data(80:end,:,:);
% % fs = 20e6;
% % speed = 1540; %m/s in body
% % pixel_size_through_depth = 0.5*(speed/fs); 
% % 
% % for ii = 1:max(size(data))
% %     time_array_all(ii) = ii/fs;
% % end
% % 
% % for cc = 1:128
% % for bb = 1:128
% %     time_array(:,bb,cc) = time_array_all;
% % end
% % end
% % 
% % channel = [[-63.5:1:63.5]];
% % 
% % for beam = 1:128
% %     
% % for jj = 1:max(size(data)) %jj=row
% %     
% % depth = jj*pixel_size_through_depth; %m
% % 
% % data_matrix = data;
% % [rows_data_matrix col_data_matrix z_data_matrix] = size(data_matrix);
% % 
% % for ii = 1:(length(channel))
% %     xe(ii) = 0.1953e-3*abs(channel(ii)); 
% %     d(ii) = (xe(ii)^2+depth^2)^0.5 + depth;
% %     time_to_point(ii) = d(ii)/speed;
% % end
% % 
% % delay_matrix(jj,:,beam) = time_to_point; %delays
% % 
% % end
% % 
% % for aa = 1:128
% %     delayed_channel(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix(1:rows_data_matrix,aa,beam),'linear');
% % end
% % 
% % 
% % end
% % 
% % 
% % for ll = 1:numel(delayed_channel)
% %     if isnan(delayed_channel(ll))==1
% %         delayed_channel(ll) = 0;
% %     end
% % end
% % 
% % figure;
% % min_data = min(min(min(delayed_channel)));
% % max_data = max(max(max(delayed_channel)));
% % imagesc(delayed_channel(:,:),[min_data, max_data])
% % colormap('gray');
% % title('Channel data with delays (anecoicCystData.mat), Problem 4');
% % 
% % summed_channels = sum(delayed_channel,2);
% % figure;
% % imagesc(20*log10(abs(hilbert(summed_channels(:,:)))));
% % colormap('gray');
% % title('Compressed B-mode image (anecoicCystData.mat), Problem 4');

clear all;
load('pointTargetData.mat');

data = veraStrct.data(80:end,:,:);
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = 0.5*(speed/fs); 

for ii = 1:max(size(data))
    time_array_all(ii) = ii/fs;
end

for cc = 1:128
for bb = 1:128
    time_array(:,bb,cc) = time_array_all;
end
end

channel = [[-63.5:1:63.5]];

%4 beams
data = veraStrct.data(80:end,:,1:4:128);
[rows_d col_d z_d] = size(data);

for beam = 1:z_d
    
for jj = 1:max(size(data)) %jj=row
    
depth = jj*pixel_size_through_depth; %m

data_matrix = data;
[rows_data_matrix col_data_matrix z_data_matrix] = size(data_matrix);

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    lat(ii) =  0.1953e-3*channel(ii); 
    d_1(ii) = ((xe(ii)-((-3/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_1(ii) = d_1(ii)/speed;
    d_2(ii) = ((xe(ii)-((-1/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_2(ii) = d_2(ii)/speed;
    d_3(ii) = ((xe(ii)-((1/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_3(ii) = d_3(ii)/speed;
    d_4(ii) = ((xe(ii)-((3/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_4(ii) = d_4(ii)/speed;
end

delay_matrix_1(jj,:,beam) = time_to_point_1; %delays
delay_matrix_2(jj,:,beam) = time_to_point_2;
delay_matrix_3(jj,:,beam) = time_to_point_3; %delays
delay_matrix_4(jj,:,beam) = time_to_point_4;

end

for aa = 1:128
    delayed_channel_1(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_1(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_2(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_2(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_3(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_3(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_4(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_4(1:rows_data_matrix,aa,beam),'linear');
end


end
ax = [1:rows_data_matrix]*pixel_size_through_depth;

[delay_rows delay_col delay_beams] = size(delayed_channel_1);
delayed_channel = zeros(rows_data_matrix,128,128);
first = 1;
second = 2;
third = 3;
fourth = 4;
for ind = 1:delay_beams %odd...delayed_channel_1
    delayed_channel(:,:,first) = delayed_channel_1(:,:,ind);
    first = first+4;
    delayed_channel(:,:,second) = delayed_channel_2(:,:,ind);
    second = second+4;
    delayed_channel(:,:,third) = delayed_channel_3(:,:,ind);
    third = third+4;
    delayed_channel(:,:,fourth) = delayed_channel_4(:,:,ind);
    fourth = fourth+4;
end

for ll = 1:numel(delayed_channel)
    if isnan(delayed_channel(ll))==1
        delayed_channel(ll) = 0;
    end
end

summed_channels = sum(delayed_channel,2);
figure;
imagesc(lat,ax,20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
axis image;
title('4 parallel beams, (pointTargetData.mat)');

clear all;
load('pointTargetData.mat');

data = veraStrct.data(80:end,:,:);
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = 0.5*(speed/fs); 

for ii = 1:max(size(data))
    time_array_all(ii) = ii/fs;
end

for cc = 1:128
for bb = 1:128
    time_array(:,bb,cc) = time_array_all;
end
end

channel = [[-63.5:1:63.5]];

%8 beams
data = veraStrct.data(80:end,:,1:8:128);
[rows_d col_d z_d] = size(data);

for beam = 1:z_d
    
for jj = 1:max(size(data)) %jj=row
    
depth = jj*pixel_size_through_depth; %m

data_matrix = data;
[rows_data_matrix col_data_matrix z_data_matrix] = size(data_matrix);

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    lat(ii) =  0.1953e-3*channel(ii); 
    d_1(ii) = ((xe(ii)-((-7/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_1(ii) = d_1(ii)/speed;
    d_2(ii) = ((xe(ii)-((-5/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_2(ii) = d_2(ii)/speed;
    d_3(ii) = ((xe(ii)-((-3/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_3(ii) = d_3(ii)/speed;
    d_4(ii) = ((xe(ii)-((-1/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_4(ii) = d_4(ii)/speed;
    d_5(ii) = ((xe(ii)-((1/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_5(ii) = d_5(ii)/speed;
    d_6(ii) = ((xe(ii)-((3/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_6(ii) = d_6(ii)/speed;
    d_7(ii) = ((xe(ii)-((5/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_7(ii) = d_7(ii)/speed;
    d_8(ii) = ((xe(ii)-((7/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_8(ii) = d_8(ii)/speed;
end

delay_matrix_1(jj,:,beam) = time_to_point_1; %delays
delay_matrix_2(jj,:,beam) = time_to_point_2;
delay_matrix_3(jj,:,beam) = time_to_point_3; %delays
delay_matrix_4(jj,:,beam) = time_to_point_4;
delay_matrix_5(jj,:,beam) = time_to_point_5;
delay_matrix_6(jj,:,beam) = time_to_point_6;
delay_matrix_7(jj,:,beam) = time_to_point_7;
delay_matrix_8(jj,:,beam) = time_to_point_8;


end

for aa = 1:128
    delayed_channel_1(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_1(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_2(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_2(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_3(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_3(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_4(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_4(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_5(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_5(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_6(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_6(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_7(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_7(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_8(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_8(1:rows_data_matrix,aa,beam),'linear');
end


end
ax = [1:rows_data_matrix]*pixel_size_through_depth;

[delay_rows delay_col delay_beams] = size(delayed_channel_1);
delayed_channel = zeros(rows_data_matrix,128,128);
first = 1;
second = 2;
third = 3;
fourth = 4;
fifth = 5;
sixth = 6;
seventh = 7;
eighth = 8;
for ind = 1:delay_beams %odd...delayed_channel_1
    delayed_channel(:,:,first) = delayed_channel_1(:,:,ind);
    first = first+8;
    delayed_channel(:,:,second) = delayed_channel_2(:,:,ind);
    second = second+8;
    delayed_channel(:,:,third) = delayed_channel_3(:,:,ind);
    third = third+8;
    delayed_channel(:,:,fourth) = delayed_channel_4(:,:,ind);
    fourth = fourth+8;
    delayed_channel(:,:,fifth) = delayed_channel_5(:,:,ind);
    fifth = fifth +8;
     delayed_channel(:,:,sixth) = delayed_channel_6(:,:,ind);
    sixth = sixth+8;
     delayed_channel(:,:,seventh) = delayed_channel_7(:,:,ind);
    seventh = seventh +8;
     delayed_channel(:,:,eighth) = delayed_channel_8(:,:,ind);
    eighth = eighth+8;
end

for ll = 1:numel(delayed_channel)
    if isnan(delayed_channel(ll))==1
        delayed_channel(ll) = 0;
    end
end

summed_channels = sum(delayed_channel,2);
figure;
imagesc(lat,ax,20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
axis image;
title('8 parallel beams, (pointTargetData.mat)');

clear all;
load('pointTargetData.mat');

data = veraStrct.data(80:end,:,:);
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = 0.5*(speed/fs); 

for ii = 1:max(size(data))
    time_array_all(ii) = ii/fs;
end

for cc = 1:128
for bb = 1:128
    time_array(:,bb,cc) = time_array_all;
end
end

channel = [[-63.5:1:63.5]];

%16 beams
data = veraStrct.data(80:end,:,1:16:128);
[rows_d col_d z_d] = size(data);

for beam = 1:z_d
    
for jj = 1:max(size(data)) %jj=row
    
depth = jj*pixel_size_through_depth; %m

data_matrix = data;
[rows_data_matrix col_data_matrix z_data_matrix] = size(data_matrix);

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    lat(ii) = 0.1953e-3*channel(ii); 
    d_1(ii) = ((xe(ii)-((-15/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_1(ii) = d_1(ii)/speed;
    d_2(ii) = ((xe(ii)-((-13/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_2(ii) = d_2(ii)/speed;
    d_3(ii) = ((xe(ii)-((-11/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_3(ii) = d_3(ii)/speed;
    d_4(ii) = ((xe(ii)-((-9/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_4(ii) = d_4(ii)/speed;
    d_5(ii) = ((xe(ii)-((-7/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_5(ii) = d_5(ii)/speed;
    d_6(ii) = ((xe(ii)-((-5/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_6(ii) = d_6(ii)/speed;
    d_7(ii) = ((xe(ii)-((-3/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_7(ii) = d_7(ii)/speed;
    d_8(ii) = ((xe(ii)-((-1/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_8(ii) = d_8(ii)/speed;
    d_9(ii) = ((xe(ii)-((1/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_9(ii) = d_9(ii)/speed;
    d_10(ii) = ((xe(ii)-((3/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_10(ii) = d_10(ii)/speed;
    d_11(ii) = ((xe(ii)-((5/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_11(ii) = d_11(ii)/speed;
    d_12(ii) = ((xe(ii)-((7/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_12(ii) = d_12(ii)/speed;
    d_13(ii) = ((xe(ii)-((9/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_13(ii) = d_13(ii)/speed;
    d_14(ii) = ((xe(ii)-((11/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_14(ii) = d_14(ii)/speed;
    d_15(ii) = ((xe(ii)-((13/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_15(ii) = d_15(ii)/speed;
    d_16(ii) = ((xe(ii)-((15/2)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_16(ii) = d_16(ii)/speed;
end

delay_matrix_1(jj,:,beam) = time_to_point_1; %delays
delay_matrix_2(jj,:,beam) = time_to_point_2;
delay_matrix_3(jj,:,beam) = time_to_point_3; %delays
delay_matrix_4(jj,:,beam) = time_to_point_4;
delay_matrix_5(jj,:,beam) = time_to_point_5;
delay_matrix_6(jj,:,beam) = time_to_point_6;
delay_matrix_7(jj,:,beam) = time_to_point_7;
delay_matrix_8(jj,:,beam) = time_to_point_8;
delay_matrix_9(jj,:,beam) = time_to_point_9;
delay_matrix_10(jj,:,beam) = time_to_point_10;
delay_matrix_11(jj,:,beam) = time_to_point_11;
delay_matrix_12(jj,:,beam) = time_to_point_12;
delay_matrix_13(jj,:,beam) = time_to_point_13;
delay_matrix_14(jj,:,beam) = time_to_point_14;
delay_matrix_15(jj,:,beam) = time_to_point_15;
delay_matrix_16(jj,:,beam) = time_to_point_16;

end

for aa = 1:128
    delayed_channel_1(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_1(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_2(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_2(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_3(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_3(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_4(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_4(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_5(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_5(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_6(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_6(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_7(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_7(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_8(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_8(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_9(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_9(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_10(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_10(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_11(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_11(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_12(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_12(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_13(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_13(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_14(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_14(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_15(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_15(1:rows_data_matrix,aa,beam),'linear');
    delayed_channel_16(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix_16(1:rows_data_matrix,aa,beam),'linear');
end


end
ax = [1:rows_data_matrix]*pixel_size_through_depth;

[delay_rows delay_col delay_beams] = size(delayed_channel_1);
delayed_channel = zeros(rows_data_matrix,128,128);
first = 1;
second = 2;
third = 3;
fourth = 4;
fifth = 5;
sixth = 6;
seventh = 7;
eighth = 8;
ninth = 9;
tenth = 10;
eleven = 11;
twelve = 12;
thirteen = 13;
fourteen = 14;
fifteen = 15;
sixteen = 16;
for ind = 1:delay_beams %odd...delayed_channel_1
    delayed_channel(:,:,first) = delayed_channel_1(:,:,ind);
    first = first+16;
    delayed_channel(:,:,second) = delayed_channel_2(:,:,ind);
    second = second+16;
    delayed_channel(:,:,third) = delayed_channel_3(:,:,ind);
    third = third+16;
    delayed_channel(:,:,fourth) = delayed_channel_4(:,:,ind);
    fourth = fourth+16;
    delayed_channel(:,:,fifth) = delayed_channel_5(:,:,ind);
    fifth = fifth +16;
     delayed_channel(:,:,sixth) = delayed_channel_6(:,:,ind);
    sixth = sixth+16;
     delayed_channel(:,:,seventh) = delayed_channel_7(:,:,ind);
    seventh = seventh +16;
     delayed_channel(:,:,eighth) = delayed_channel_8(:,:,ind);
    eighth = eighth+16;
    delayed_channel(:,:,ninth) = delayed_channel_9(:,:,ind);
    ninth = ninth+16;
    delayed_channel(:,:,tenth) = delayed_channel_10(:,:,ind);
    tenth = tenth+16;
    delayed_channel(:,:,eleven) = delayed_channel_11(:,:,ind);
    eleven = eleven+16;
    delayed_channel(:,:,twelve) = delayed_channel_12(:,:,ind);
    twelve = twelve+16;
    delayed_channel(:,:,thirteen) = delayed_channel_13(:,:,ind);
    thirteen = thirteen +16;
     delayed_channel(:,:,fourteen) = delayed_channel_14(:,:,ind);
    fourteen = fourteen+16;
     delayed_channel(:,:,fifteen) = delayed_channel_15(:,:,ind);
    fifteen = fifteen +16;
     delayed_channel(:,:,sixteen) = delayed_channel_16(:,:,ind);
    sixteen = sixteen+16;
end

for ll = 1:numel(delayed_channel)
    if isnan(delayed_channel(ll))==1
        delayed_channel(ll) = 0;
    end
end

summed_channels = sum(delayed_channel,2);
figure;
imagesc(lat,ax,20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
axis image;
title('16 parallel beams, (pointTargetData.mat)');

