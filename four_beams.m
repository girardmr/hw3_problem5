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


%4 beams
data = veraStrct.data(80:end,:,1:4:128);
[rows_d col_d z_d] = size(data);

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

for beam = 1:z_d
    
for jj = 1:max(size(data)) %jj=row
    
depth = jj*pixel_size_through_depth; %m

data_matrix = data;
[rows_data_matrix col_data_matrix z_data_matrix] = size(data_matrix);

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    d_1(ii) = ((xe(ii)-((-3/4)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_1(ii) = d_1(ii)/speed;
    d_2(ii) = ((xe(ii)-((-1/4)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_2(ii) = d_2(ii)/speed;
    d_3(ii) = ((xe(ii)-((1/4)*0.1953e-3))^2+depth^2)^0.5 + depth;
    time_to_point_3(ii) = d_3(ii)/speed;
    d_4(ii) = ((xe(ii)-((3/4)*0.1953e-3))^2+depth^2)^0.5 + depth;
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
imagesc(20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
title('4 parallel beams, (pointTargetData.mat)');


