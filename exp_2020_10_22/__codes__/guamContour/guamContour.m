close all;
name = "Triton Contour"; %Asan Beach
% info = shapeinfo('guam-topography-contours-10m/contours_10m.shp');
% S = shaperead('guam-topography-contours-10m/contours_10m.shp');
% mapshow(S,'LineWidth',2.5);
% hold on;

%Find the data file here
%https://pae-paha.pacioos.hawaii.edu/thredds/ncss/usgs_dem_10m_guam?var=elev&horizStride=1
ncfile = 'usgs_dem_10m_guam.nc';
ncinfo(ncfile);
ncdisp(ncfile);

 elevation = ncread(ncfile,'elev');
 lon = ncread(ncfile,'lon');
 lat = ncread(ncfile,'lat');

 lonGrid = repmat(lon,1,size(elevation,2));
 latGrid = transpose(lat);
 latGrid = repmat(latGrid,size(elevation,1),1);
 
 %To bound the lon and lat to make the code run faster
 boundLon = lon(4074-50:4074+50);
 boundLat = lat(3481-50:3481+50);
 boundElev = elevation(4074-50:4074+50,3481-50:3481+50);
 
 %Lat and Lon need to be the same size as elev for contourm to work
 boundLon = repmat(boundLon,1,size(boundElev,2));
 boundLat = transpose(boundLat);
 boundLat = repmat(boundLat,size(boundElev,1),1);
%boundElev = round(boundElev);

%[10,20,30,40],
[C,h] = contourm(boundLat,boundLon,boundElev,[125,130,135,140,145,150],'k');
ctext = clabelm(C,h,'manual');
set(ctext,'FontSize',18,'Margin',0.1); 
hold on;
%contourf(elevation);

 ylim([13.532259-0.004 13.532259+0.004])   
 xlim([144.872159-0.003 144.872159+0.004])



% [X,Y] = ll2utm(13.473615,144.709086); %Nat Park Asan Beach 2847 2314
%[X,Y] = ll2utm(13.532259,144.872159); %University of Guam Triton Farm 3481 4074
% xlim([X-500 X+500])   %Set map to 1km box around center of site
% ylim([Y-500 Y+500])
 
% Asan Beach CRB Traps
% x = [13.473947 13.473172 13.475753 13.474981 13.474682 13.473682 13.473717 13.471296];
% y = [144.708133 144.708293 144.707848 144.708669 144.708796 144.711611 144.712720 144.706369];
 
%Triton Farm CRB Traps
 x = [13.532153 13.530453 13.531464 13.529722 13.530515 13.529263 13.532215 13.533227 13.535053];
 y = [144.873579 144.872495 144.871648 144.871551 144.870457 144.870446 144.870403 144.873718 144.873901];
 
%Convert to Northing/Easting
%  [x,y] = ll2utm(x,y);
 plot(y,x,'k.','MarkerSize',28);
title('Asan Beach Contour Map');
xlabel('Lon (°)',"FontSize",12);
ylabel('Lat (°)',"FontSize",12);

% Assan beach arrows                         
% p1 = [2.5189e5, 1.49035e6];         % First Point 
% p2 = [2.5169e5, 1.4902e6];          % Second Point
% p3 = [2.51934e5,1.4902e6];
% dp = p2-p1;                         % Difference
% dp1 = p3-p1;
% quiver(p1(1),p1(2),dp(1),dp(2),0,'linewidth',5,'color',[1,0,0]);
% quiver(p1(1),p1(2),dp1(1),dp1(2),0,'linewidth',5,'color',[1,0,0]);

legend('CRB Trap Locations','color','white','Location','southeast','FontSize', 24)

% Text used for Assan beach
% txt = " Max height: 40m";
 contourtxt = "Contour lines: 5m";
%  ctxt = text(144.711,13.477,contourtxt,"FontSize",24);
 ctxt = text(144.8744,13.5291,contourtxt,"FontSize",24);

 %text(2.521e5,1.491075e6,contourtxt,"FontSize",24);
 
helperAdjustFigure(gcf,name);



%ylim([242102.465087609;278942.845662924])   %Nat Park Asan Beach
%xlim([1464222.63864937;1510526.99747624])

% xlim([144.872159-0.025 144.872159+0.025])   %University of Guam Triton Farm
% ylim([13.532259-0.025 13.532259+0.025])