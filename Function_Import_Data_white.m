function f = Function_Import_Data_white(file)

max_Exp_time=19.25;     %24.5

filename = file;
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
generation = dataArray{:, 1};
cell0 = dataArray{:, 2};
birthwhitestart = dataArray{:, 3};
whiteendgreenstart = dataArray{:, 4};
greenend = dataArray{:, 5};
roundstartMstart = dataArray{:, 6};
celldivision = dataArray{:, 7};
celldeath = dataArray{:, 8};
celllost = dataArray{:, 9};
survived = dataArray{:, 10};
time = dataArray{:, 11};

%%%%% Exp_end: cell survived exeriment
survived((celldeath-time)*0.25>max_Exp_time)=1;
celldeath((celldeath-time)*0.25>max_Exp_time)=nan;

clearvars filename delimiter startRow formatSpec fileID dataArray ans;


t_Birth_TRAIL_white=[];       %TRAIL after Birth
t_TRAIL_Death_white=[];       %TRAIL to death
t_TRAIL_WhiteEnd_F0_white=[]; %after TRAIL
t_TRAIL_Div_F0_white=[];      %after TRAIL
t_TRAIL_WhiteEnd_F1_white=[]; %after TRAIL
numb=[];

for i=1:length(generation)
    if generation(i) == 0 && ~isnan(celldeath(i)) && ~isnan(birthwhitestart(i))   %cells died in F0, birth is known
        if whiteendgreenstart(i) < celldeath(i) && ~isnan(whiteendgreenstart(i))   % cells died in SG2M
            t_Birth_TRAIL_white=[t_Birth_TRAIL_white, (time(i)-birthwhitestart(i))*0.25];
            t_TRAIL_Death_white=[t_TRAIL_Death_white, (celldeath(i)-time(i))*0.25];
            t_TRAIL_WhiteEnd_F0_white=[t_TRAIL_WhiteEnd_F0_white, (whiteendgreenstart(i)-time(i))*0.25];
            t_TRAIL_Div_F0_white=[t_TRAIL_Div_F0_white, nan];
            t_TRAIL_WhiteEnd_F1_white=[t_TRAIL_WhiteEnd_F1_white, nan];
            numb=[numb,cell0(i)];
        else %cells died in G1
            t_Birth_TRAIL_white=[t_Birth_TRAIL_white, (time(i)-birthwhitestart(i))*0.25];
            t_TRAIL_Death_white=[t_TRAIL_Death_white, (celldeath(i)-time(i))*0.25];
            t_TRAIL_WhiteEnd_F0_white=[t_TRAIL_WhiteEnd_F0_white, nan];
            t_TRAIL_Div_F0_white=[t_TRAIL_Div_F0_white, nan];
            t_TRAIL_WhiteEnd_F1_white=[t_TRAIL_WhiteEnd_F1_white, nan];  
            numb=[numb,cell0(i)];
        end
    elseif generation(i) == 1 && ~isnan(celldeath(i)) %cells died in F1
        for m=1:4
                if round(cell0(i))==round(cell0(i-m)) && generation(i-m)==0      %i-m was mother
                    if whiteendgreenstart(i) < celldeath(i) && ~isnan(whiteendgreenstart(i))   % cells died in SG2M F1
                        t_Birth_TRAIL_white=[t_Birth_TRAIL_white, (time(i)-birthwhitestart(i-m))*0.25];  
                        t_TRAIL_Death_white=[t_TRAIL_Death_white, (celldeath(i)-time(i))*0.25];
                        t_TRAIL_WhiteEnd_F0_white=[t_TRAIL_WhiteEnd_F0_white, (whiteendgreenstart(i-m)-time(i))*0.25];
                        t_TRAIL_Div_F0_white=[t_TRAIL_Div_F0_white, (celldivision(i-m)-time(i))*0.25];
                        t_TRAIL_WhiteEnd_F1_white=[t_TRAIL_WhiteEnd_F1_white, (whiteendgreenstart(i)-time(i))*0.25];    
                        numb=[numb,cell0(i)];
                    else    %cells died in G1 F1
                        t_Birth_TRAIL_white=[t_Birth_TRAIL_white, (time(i)-birthwhitestart(i-m))*0.25];  
                        t_TRAIL_Death_white=[t_TRAIL_Death_white, (celldeath(i)-time(i))*0.25];
                        t_TRAIL_WhiteEnd_F0_white=[t_TRAIL_WhiteEnd_F0_white, (whiteendgreenstart(i-m)-time(i))*0.25];
                        t_TRAIL_Div_F0_white=[t_TRAIL_Div_F0_white, (celldivision(i-m)-time(i))*0.25];
                        t_TRAIL_WhiteEnd_F1_white=[t_TRAIL_WhiteEnd_F1_white, nan];
                        numb=[numb,cell0(i)];
                    end
                end
        end
    end
end

        
   
%%%%% exclude cells for which beginning of phase is not known %%%%%%%%%%%%%
f= [t_Birth_TRAIL_white(~isnan(t_Birth_TRAIL_white));t_TRAIL_WhiteEnd_F0_white(~isnan(t_Birth_TRAIL_white));t_TRAIL_Div_F0_white(~isnan(t_Birth_TRAIL_white));t_TRAIL_WhiteEnd_F1_white(~isnan(t_Birth_TRAIL_white));t_TRAIL_Death_white(~isnan(t_Birth_TRAIL_white))];  

end



    
    
