load('../../data/CRB/CRB_Traps.mat', 'CRB_Traps');
trap_coordinates = [CRB_Traps.LATITUDE, CRB_Traps.LONGITUDE]
cellLats = [];
cellLons = [];
Trap_WAIKIKI = [];
Wifi_TWO = [];
Wifi_THREE = [];
Wifi_FOUR = [];

out = inpolygon(trap_coordinates(:,1),trap_coordinates(:,2),[21.27305;21.271125;21.275695;21.278473;21.280912;21.278864;21.282558;21.287785;21.289175;],[-157.817059;-157.824017;-157.827495;-157.835877;-157.838163;-157.842426;-157.851271;-157.847340;-157.835701]);
size(out)
for j=1:size(out)
    if out(j) == 1
        Trap_WAIKIKI = [Trap_WAIKIKI;trap_coordinates(j,:)];
    end
end
% out = inpolygon(trap_coordinates(:,1),trap_coordinates(:,2),[21.371434;21.302260;21.3072],[-158.026;-158.053837;-158.00]);
% for j=1:size(out)
%     if out(j) == 1
%         Wifi_TWO = [Wifi_TWO;trap_coordinates(j,:)];
%     end
% end
% out = inpolygon(trap_coordinates(:,1),trap_coordinates(:,2),[21.336;21.3712;21.3712;21.3362],[-157.958;-157.958;-157.9205;-157.9279]);
% for j=1:size(out)
%     if out(j) == 1
%         Wifi_THREE = [Wifi_THREE;trap_coordinates(j,:)];
%     end
% end 
% out = inpolygon(trap_coordinates(:,1),trap_coordinates(:,2),[21.3433;21.33554;21.325;21.325;21.3376],[-157.92175;-157.926;-157.925;-157.913;-157.912]);
% for j=1:size(out)
%     if out(j) == 1
%         Wifi_FOUR = [Wifi_FOUR;trap_coordinates(j,:)];
%     end
% end  
size(Trap_WAIKIKI)
Trap_WAIKIKI
% size(Wifi_TWO)
% size(Wifi_THREE)
% size(Wifi_FOUR)
save('Trap_WAIKIKI.mat', 'Trap_WAIKIKI');
% save('Wifi_TWO.mat', 'Wifi_TWO');
% save('Wifi_THREE.mat', 'Wifi_THREE');
% save('Wifi_FOUR.mat', 'Wifi_FOUR');