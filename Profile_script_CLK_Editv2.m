%% 
%
%
%
%
% Requirements MATLAB 2016a+ Mapping Toolbox 
%
%
%

% General Comments 
% Assumes spehrical Mars ellipspod with radius 3,396,190 m 
% https://astrogeology.usgs.gov/search/map/Mars/Topography/HRSC_MOLA_Blend/Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2
% 
%

%% settings 
DEMfile_geotiff='NOCTIS_HRSC_MOLA_BlendDEM_Global_200mp';  % MOLA DEM in geotiff format 
grabencenterline_shp='NL_DP_18.shp'; % center line trace of graben floor 
fno=18;  % this is the fault number to read in & used for file name 
L_km=6;  % this is the half length of the fault normal profile centered on trace
S_km=2;  % this is the spacing between fault normal profiles 

%% read in data for a shape file having only on fault trace 
S=shaperead(grabencenterline_shp);
lat=S.Y; lat=lat(1:end-1); % name S.Y lat and take all but the last value, NAN 
lon=S.X; lon=lon(1:end-1); % name S.X lon and take all but the last value, NAN 

%%  read the MOLA elevation data 
[Z, R] = geotiffread(DEMfile_geotiff); 

%% trim and plot the MOLA DEM to around the fault
tmp1=range(lat); tmp2=range(lon)/2; 
LATLIM=[min(lat)-km2deg(5+L_km,3396.190), max(lat)+km2deg(5+L_km,3396.190)]; 
LONLIM=[min(lon)-km2deg(5+L_km,3396.190), max(lon)+km2deg(5+L_km,3396.190)]; 
[ZT, RT] = maptrims(Z,R,LATLIM,LONLIM); ZT=double(ZT);
clims=quantile(ZT(:),[0.25,0.75]);  % find quatiile elevations in AOI  
clear Z R;

%% Resample graben trace at S_km spacing 
% densifies points along trace to S_km spacing 
sp=km2deg(S_km,3396.190); [ptlat,ptlon]=interpm(lat,lon,sp); 

% calculate the cumulative distance along fault in km
        tmpd=size(ptlat); 
        for d=1:length(ptlat)-1 
            tmpd(d)=distance(ptlat(d),ptlon(d),ptlat(d+1),ptlon(d+1),[3396.190 0]);
        end
        tmpd=cat(2,0,tmpd); km_along_fault=cumsum(tmpd);

 %calculate the cumulative distance along fault using original points
         tmpd=nan(1,length(lat)-1); 
         for d=1:length(lat)-1 
             tmpd(d)=distance(lat(d),lon(d),lat(d+1),lon(d+1),[3396.190 0]);
         end
         tmpd=cat(2,0,tmpd); km_along_fault=cumsum(tmpd);
     
 % resample at exactly S_km based on distance 
         ptsS=0:S_km:km_along_fault(end); % points to sample (plus end point)
         ptlat=interp1(km_along_fault,lat,ptsS); ptlat=cat(2,ptlat,lat(end));
         ptlon=interp1(km_along_fault,lon,ptsS); ptlon=cat(2,ptlon,lon(end));
         
 % now recalculate distance for new points at S_km spacing 
         tmpd=nan(1,length(ptlat)-1); 
         for d=1:length(ptlat)-1 
             tmpd(d)=distance(ptlat(d),ptlon(d),ptlat(d+1),ptlon(d+1),[3396.190 0]);
         end
         tmpd=cat(2,0,tmpd); km_along_fault=cumsum(tmpd);

%% azimuth of the fault at every point 
% azimuth between points 
    AZ=nan(size(ptlat)); 
    for a=2:length(ptlat)-1  
        AZ(a)=azimuth(ptlat(a-1),ptlon(a-1),ptlat(a+1),ptlon(a+1),[3396.190 0],'degrees'); 
    end
% smooth the azimuth data 
    smfac=floor(length(ptlat)/6); if rem(smfac,2)==0; smfac=smfac+1; end   % make sure its an odd number 
    AZ(1)=AZ(2);  AZ(end)=AZ(end-1);  AZ=smoothazimuth(AZ,smfac); 

    
% %% Now extract and plot for selection 
dataQuest=input('Did you save the previous datafile out? y or n ','s');
 if strcmp(dataQuest,'y')
     dataout=nan(length(ptlat),22);  
 else
     disp('Please save data out and rerun code')
    return
 end


 %% MAP and Profile 
figure('Position', [100, 100, 1450, 800]) 
i=1; % start with the first profile 
while i <= length(ptlat)  
    
   subplot(2,1,1); hold off; axesm('MapProjection','Robinson','Geoid',...
       [3396190 0]); h0=meshm(ZT, RT); caxis(clims);   % colored bathy 
   hold on; h1=plotm(lat,lon,'k','LineWidth',0.5);     %original line 
   h2=plotm(ptlat(i),ptlon(i),'.m','MarkerSize',12); hold on; % points at S_km spacing 
   
   % trace of fault & current profile 
    h3=plotm(dataout(:,8),dataout(:,9),'w','LineWidth',1);  h4=plotm(dataout(:,11),dataout(:,12),'w','LineWidth',1); 
    title(['Fault #: ' num2str(fno) '  Profile #' num2str(i) '']); 
 
   % calculate the profile and display on map 
     tmp=track1([ptlat(i),ptlat(i)]',[ptlon(i),ptlon(i)]',[wrapTo360(AZ(i)+90),wrapTo360(AZ(i)-90)]',[L_km,L_km]',[3396.190 0],'degrees',1); % find end points 
     tmp=sortrows([tmp(1),tmp(3); tmp(2), tmp(4)],1); % sort them so that the most southern point is first 
     [R_z,R_r,R_lat,R_lon] = mapprofile(ZT,RT,tmp(:,1),tmp(:,2),[3396.190 0],'gc','nearest');   % this will extract the profile between the end pints 
     T = table(R_r, R_z, R_lat, R_lon) ;
     
   % write out topoprofile to .csv file
     fno_s=sprintf('%03.0f',fno) ;
     % name .csv according to fault number and profile numer
     filename = strcat('faultNum', fno_s, '_profile', num2str(i) , '.csv') ;
     writetable(T, filename) ;
     
     % plot the profile on the map 
     subplot(2,1,1); hold on; h5=plotm(R_lat,R_lon,'k','LineWidth',1); h6=plotm(R_lat(1),R_lon(1),'ro','MarkerSize',8);  
 
     % plot the range verses elevation data for the profile and add the local slope 
     subplot(2,1,2); plot(R_r,R_z); hold on; plot(R_r(1),R_z(1),'ro','MarkerSize',12); % profile data 
     scatter(R_r(2:end),R_z(2:end),8,smooth(atand(diff(R_z/1000)./(diff(R_r))),3)); c=colorbar;  c.Label.String='Slope Deg'; c.Label.FontSize=18;    % slope 
     plot(ones(1,2).*range(R_r)/2, [min(R_z),max(R_z)],'m','LineWidth',1);   % draw a line where the trace was picked. 
     grid on; xlim([0,2*L_km]);  % set the plot limits. 
 
     % pick the points and plot them 
     title('CLICK first scarp base and then first scarp top; REPEAT FOR SECOND base & top; then ENTER','FontSize',15)  % give instructions 
     datain=ginput; %  get the point data 
     plot(datain(1,1),datain(1,2),'k+','MarkerSize',12,'LineWidth',2); text(datain(1,1)+0.3,datain(1,2),'B1','FontSize',18); % plot first base point on profile 
     plot(datain(2,1),datain(2,2),'k+','MarkerSize',12,'LineWidth',2); text(datain(2,1)+0.3,datain(2,2),'T1','FontSize',18); % plot first top point on profile 
     plot(datain(3,1),datain(3,2),'k+','MarkerSize',12,'LineWidth',2); text(datain(3,1)+0.3,datain(3,2),'B2','FontSize',18); % plot second base point on profile
     plot(datain(4,1),datain(4,2),'k+','MarkerSize',12,'LineWidth',2); text(datain(4,1)+0.3,datain(4,2),'T2','FontSize',18); hold off % plot second top point on profile
     mkB1=find(R_r < datain(1,1),1,'last');   mkT1=find(R_r > datain(2,1),1,'first');  mkB2=find(R_r < datain(3,1),1,'last'); mkT2=find(R_r > datain(4,1),1,'first'); % get there lat and lons 
     subplot(2,1,1); h7=plotm(R_lat(mkB1),R_lon(mkB1),'w.','MarkerSize',12); h8=plotm(R_lat(mkT1),R_lon(mkT1),'w.','MarkerSize',12); % plot them on the map (1/2)
     subplot(2,1,1); h9=plotm(R_lat(mkB2),R_lon(mkB2),'w.','MarkerSize',12); h10=plotm(R_lat(mkT2),R_lon(mkT2),'w.','MarkerSize',12); % plot them on map map (2/2)
     
     clims=[min(R_z),max(R_z)];
     % now confirm you choices 
     subplot(2,1,2);
     title('Enter 1 to Accept, 2 to repeat, or  0 to skip this profile','FontSize',20)
     resp=input('Enter 1 to Accept, 2 to repeat, or  0 to skip this profile:    '); 

  % now do different things depending on if you want accept, repeat or skip
  % the profile 
    if resp==1
       disp('saving this profile') 
       throw1=abs(R_z(mkT1)-R_z(mkB1));
       throw2=abs(R_z(mkT2)-R_z(mkB2));
       baseWidth=distance(R_lat(mkB2), R_lon(mkB2), R_lat(mkB1), R_lon(mkB1),[3396190 0]);
       strainWidth=distance(R_lat(mkT2), R_lon(mkT2), R_lat(mkT1), R_lon(mkT1),[3396190 0]);
       strainAmount=(baseWidth/strainWidth)*100;
       slopeangle1=atand((R_z(mkT1)-R_z(mkB1))/distance(R_lat(mkB1), R_lon(mkB1), R_lat(mkT1), R_lon(mkT1),[3396190 0]));  % inverse tan of elevation in m over distance in meters 
       slopeangle2=atand((R_z(mkT2)-R_z(mkB2))/distance(R_lat(mkB2), R_lon(mkB2), R_lat(mkT2), R_lon(mkT2),[3396190 0]));
       tmpout=[i, km_along_fault(i), AZ(i), throw1, throw2, slopeangle1, slopeangle2, R_lat(mkB1), R_lon(mkB1), R_z(mkB1), R_lat(mkT1), R_lon(mkT1), R_z(mkT1), R_lat(mkB2), R_lon(mkB2), R_z(mkB2), R_lat(mkT2), R_lon(mkT2), R_z(mkT2), baseWidth, strainWidth, strainAmount]; 
       dataout(i,:)=tmpout; 
       delete([h5,h6,h7,h8,h9,h10])
       i=i+1; % can use -1 if you want to go the opposite way along the fault
    elseif resp==0
        i=i+1; % can use -1 if you want to go the opposite way along the fault
        disp('not recoding a throw for this profile') 
        delete([h5,h6,h7,h8,h9,h10])
    else 
        disp('repeating the profile') 
        delete([h5,h6,h7,h8,h9,h10])
    end
end


%% now write out the data into a text file 
    fno_s=sprintf('%03.0f',fno);
    % writes file with format fault_no***_v*dat/time inserted here*
    fileout=strcat('faultNum_', fno_s, '_v', datestr(now,'ddmmyyyy'), '_', datestr(now,'HHMM'), '.csv'); 
    fid=fopen(fileout,'w'); 
    fprintf(fid, 'profileNo, km_along_fault, fault_azimuth, throw_m, throw_m2, slope_deg, slope_deg2, lat_base1, lon_base1, elev_base1, lat_top1, lon_top1, elev_top1, lat_base2, lon_base2, elev_base2, lat_top2, lon_top2, elev_top2, baseWidth_km, strainWidth_km, strainAmount\n');
    fprintf(fid,'%03.0f, %1.3f, %1.3f, %1.3f, %1.1f, %15.12f, %15.12f, %1.1f, %15.12f, %15.12f, %1.1f, %1.3f, %1.1f, %15.12f, %15.12f, %1.1f, %15.12f, %15.12f, %1.1f, %1.3f, %1.3f, %2.0f\n', dataout'); 
    fclose(fid); 



%% plot the data throws vs dist to check data
    figure;
    colors = {[0.2 0.2 0.2]; [0.7 0.7 0.7]; [0.87 0.77 0.57]};
    plot(dataout(:,2),dataout(:,4), 'Color', colors{1}); hold on
    plot(dataout(:,2),dataout(:,5), '-.','Color', colors{2});
    xlabel('Distance(km)');
    ylabel('Throw (m)');
    yyaxis right
    ylabel('Graben Width (km)');
    plot(dataout(:,2),dataout(:,20),'bo', 'Color', colors{3});
    title('Throw Profile','FontSize',20);  hold off


% FUNCTIONS 
    function [smAZ]=smoothazimuth(AZ,smfac) 
    % AZ = vector of azimuths in degrees 
    % smfac smoothing factor, should be an odd integer. e.g., 3,7,9 
    % smAZ = running average vector mean aziuth. 

    azB=buffer(AZ,smfac,smfac-1);  % reshape the input into overlapping columns 
    for i=1:smfac-1
    azB(1:smfac-i,i)=AZ(i); %deal with first few columns that get zero padded  
    end

    % now call vector mean for each column and store the output. 
    smAZ=nan(size(AZ)); 
    for j=1:length(AZ)
    smAZ(j)=vectormean(azB(:,j)); 
    end
    end


   function meanA=vectormean(azdeg)
    cosPhase = sum(cosd(azdeg));
    sinPhase = sum(sind(azdeg));

        if (sinPhase >=0 && cosPhase >=0) 
            meanA = atand(sinPhase/cosPhase);  
        end
        if (cosPhase < 0) 
            meanA = atand(sinPhase/cosPhase) + 180;
        end
        if (sinPhase <0 && cosPhase >=0) 
            meanA = atand(sinPhase/cosPhase)+360;
        end
   end
   

