function [ max_sample_freq ] = Comparator_Max_Freq( vop, vom, vclock, num_clock_sims )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


for i=1:num_clock_sims
    vdiff.V=-(vop.V(:,i)-vom.V(:,i));
    vdiff.time=vop.time(:,i);
    
    search_1=0.5;
    search_2=0.9;
    Vdd=1.2;
    % find index = 0.1
    [~,start_index]=min(abs(Vdd*search_1-vclock.V(:,i)));
    [~,stop_index]=min(abs(Vdd*search_2-vdiff.V));
    min_dec_time(i)=vdiff.time(stop_index)-vdiff.time(start_index);
    max_freq_clock(i)=(1/min_dec_time(i))/2;
    max_sample_freq(i)=max_freq_clock(i)/12;
    fprintf('Min decision time: %4.2f ps\n',min_dec_time(i)*1e12);
    fprintf('Max sample frequency: %4.2f MHz\n',max_sample_freq(i)*1e-6);
    clearvars vdiff.V
end

end

