load('../../data/CRB/PEARL_HARBOR.mat', 'Pearlharbor')

cellLats = [];
cellLons = [];
Wifi_ONE = [];
Wifi_TWO = [];
Wifi_THREE = [];
Wifi_FOUR = [];

for i = 1:numCellSites
   cellLats = [cellLats; Pearlharbor(i).trilat];
   cellLons = [cellLons; Pearlharbor(i).trilong];
end

wifi_coordinates = [cellLats cellLons];
out = inpolygon(wifi_coordinates(:,1),wifi_coordinates(:,2),[21.329434;21.402260;21.4402;21.378],[-158.020828;-158.05;-158.00;-157.931]);
for j=1:size(out)
    if out(j) == 1
        Wifi_ONE = [Wifi_ONE;wifi_coordinates(j,:)];
    end
end
out = inpolygon(wifi_coordinates(:,1),wifi_coordinates(:,2),[21.371434;21.302260;21.3072],[-158.026;-158.053837;-158.00]);
for j=1:size(out)
    if out(j) == 1
        Wifi_TWO = [Wifi_TWO;wifi_coordinates(j,:)];
    end
end
out = inpolygon(wifi_coordinates(:,1),wifi_coordinates(:,2),[21.336;21.3712;21.3712;21.3362],[-157.958;-157.958;-157.9205;-157.9279]);
for j=1:size(out)
    if out(j) == 1
        Wifi_THREE = [Wifi_THREE;wifi_coordinates(j,:)];
    end
end 
out = inpolygon(wifi_coordinates(:,1),wifi_coordinates(:,2),[21.3433;21.33554;21.325;21.325;21.3376],[-157.92175;-157.926;-157.925;-157.913;-157.912]);
for j=1:size(out)
    if out(j) == 1
        Wifi_FOUR = [Wifi_FOUR;wifi_coordinates(j,:)];
    end
end  
size(Wifi_ONE)
size(Wifi_TWO)
size(Wifi_THREE)
size(Wifi_FOUR)
save('Wifi_ONE.mat', 'Wifi_ONE');
save('Wifi_TWO.mat', 'Wifi_TWO');
save('Wifi_THREE.mat', 'Wifi_THREE');
save('Wifi_FOUR.mat', 'Wifi_FOUR');