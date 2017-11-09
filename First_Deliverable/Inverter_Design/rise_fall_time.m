function [ rise_time, fall_time ] = rise_fall_time( data, Vdd )

% eliminate < 0 and > Vdd data
TF=data(:,2)>Vdd | data(:,2)<0;
data(TF,:)=[];

% switching occurs near 250ps and 500ps
[~,rise_index_1]=min(abs(data(:,1)-0e-12));
[~,rise_index_2]=min(abs(data(:,1)-249e-12));
[~,fall_index_1]=min(abs(data(:,1)-251e-12));
[~,fall_index_2]=min(abs(data(:,1)-499e-12));

% rise time define as between 10% and 90% of the inverter swing
test_1=0.1*Vdd;
test_2=0.9*Vdd;

[~,rise_index_low]=min(abs(data(rise_index_1:rise_index_2,2)-test_1));
[~,rise_index_high]=min(abs(data(rise_index_1:rise_index_2,2)-test_2));
rise_time=data(rise_index_high,1)-data(rise_index_low,1);

[~,fall_index_low]=min(abs(data(fall_index_1:fall_index_2,2)-test_1));
[~,fall_index_high]=min(abs(data(fall_index_1:fall_index_2,2)-test_2));
fall_time=data(fall_index_low,1)-data(fall_index_high,1);



end

