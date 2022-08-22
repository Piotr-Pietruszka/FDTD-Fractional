
% x_l = get(gca, 'XLim');
% % xlim([x_l(1), x_l(2)/2]); 
%% scale to get narrow
% f = gcf;
% position_org = f.Position;
% f.Position(4) = position_org(4)*0.7 % scale height - smaller
% % % f.Position(4) = position_org(4)*1.2 % scale height - greater
% 
% f.Position(3) = position_org(3)*1.8 % scale width - greater
% 



%% for horizontallly placed multiple plots (greater height)
x_l = get(gca, 'XLim');
xlim([x_l(1), x_l(2)* 1/2]); % signal at the begiining 
% xlim([x_l(2)* 1/2, x_l(2)]); % signal at the end 
ylim([-1, 1]);
f = gcf;
position_org = f.Position;
f.Position(4) = position_org(4)*1.2 % scale height - greater
f.Position(2) = 200; % get it down
%% scale to get narrow - movie
% f = gcf;
% position_org = f.Position;
% f.Position(4) = position_org(4)*0.7 % scale height - smaller
% f.Position(3) = position_org(3)*2 % scale width - greater
% f.Position(4) = position_org(4)*0.8 % scale height - smaller
% f.Position(3) = position_org(3)*1.5 % scale width - greater