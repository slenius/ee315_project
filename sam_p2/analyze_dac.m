clear all
close all

addpath('../common')

% given params
c = load_project_const();

load('dacdata_ab_2.mat');

fs = 50e6;
ts = 1 / fs;
dly = ts;
tbit = ts/12;



bit_time = [];
for i = 1:10
    bit_time(i) = (2*dly - tbit) - i*tbit;
end

bits_to_test = 1:9;



figure(1);
title('Differential Mode Voltage');
figure(2);
title('Common Mode Voltage');
leg = cell(1);

[wu, order] = sort(wu);

colors = jet(length(data));

for i = 1:length(data)
   di = order(i);
   t{i} = data{di}.dacout_2_p.time;
   yp = data{di}.dacout_2_p.V;
   ym = data{di}.dacout_2_m.V;
   dm{i} = yp - ym;
   cm{i} = (yp + ym) / 2;
   
   figure(1);
   plot(t{i}, dm{i}, '-', 'Color', colors(i,:));
   hold on;
   leg{i} = sprintf('wu = %0.2f um', wu(i)*1e6);
   
   figure(2);
   plot(t{i}, cm{i}, '-', 'Color', colors(i,:));
   hold on;
   
   % compute delay for each bit
   for j = 1:length(bits_to_test)
     b = bits_to_test(j);
     tref_pre = bit_time(b);
     tref_post = tref_pre + 0.85*tbit;
   
     % get the indicies of the pre-switch and post-switch time
     [~, xi_pre] = min(abs(t{i}-tref_pre));
     [~, xi_post] = min(abs(t{i}-tref_post));
   
     % get the reference voltage levels- pre, post and 99% settling
     v_pre = dm{i}(xi_pre);
     v_post = dm{i}(xi_post);
     v_ref = 0.9 * (v_post - v_pre) + v_pre;
   
     % find the voltage index for the crossing point within half lsb
     %[~, yi_ref] = min(abs(dm(i,:)(xi_pre:xi_post)-v_ref));
     %yi_ref = yi_ref + xi_pre;
     r = c.lsb / 2;
     for k = xi_post:-1:xi_pre
         a = abs(dm{i}(k) - v_post);
         if a >= r
             if dm{i}(k) > v_post
                 xq = v_post + r;
             else
                 xq = v_post - r;
             vi = t{i}(k:xi_post);
             end
             xi = dm{i}(k-2:k+2);
             ti = t{i}(k-2:k+2);
             tq = interp1(xi, ti, xq);
             
             
             % compute the settling time
             settle_time(j, i) = tq - t{i}(xi_pre);
             
             % annotate
             figure(1);
             s = sprintf('<- wu = %0.3f um, ts = %0.3f ns', wu(i) * 1e6, settle_time(j, i) * 1e9);
             h = text(tq,xq,s);
             set(h,'Rotation',90);
             break
         end
     end     
     
     bit(j) = b;
   end
   
end

[l,~] = size(settle_time);

figure(1);
legend(leg);
figure(2);
legend(leg);





leg = cell(1);
figure(3);
colors = jet(length(bits_to_test));
for i = 1:l
    leg{i} = sprintf('Bit %d', bit(i));
    plot(wu*1e6, settle_time(i,:)*1e9, '*-', 'Color', colors(i,:));
    hold on;
end
legend(leg);
ylim([0, 1]);
xlabel('Unit Width (um)')
ylabel('Delay (ns)');
title('Time to settle to 0.5 LSB');


figure(4);
plot(wu*1e6, settle_time(9,:)*1e9, 'k-*');
ylim([0, 1]);
xlabel('Unit Width (um)')
ylabel('Delay (ns)');
title('Time to settle to 0.5 LSB');
