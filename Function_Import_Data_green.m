function f=Function_Import_Data_green(file)

max_Exp_time=19.25;     %24.5


filename = file;
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
generation1 = dataArray{:, 1};
cell1 = dataArray{:, 2};
birthwhitestart1 = dataArray{:, 3};
whiteendgreenstart1 = dataArray{:, 4};
greenend1 = dataArray{:, 5};
roundstartMstart1 = dataArray{:, 6};
celldivision1 = dataArray{:, 7};
celldeath1 = dataArray{:, 8};
celllost1 = dataArray{:, 9};
survived1 = dataArray{:, 10};
time = dataArray{:, 11};

%%%%% Exp_end: cell survived exeriment
survived1((celldeath1-time)*0.25>max_Exp_time)=1;
celldeath1((celldeath1-time)*0.25>max_Exp_time)=nan;

clearvars filename delimiter startRow formatSpec fileID dataArray ans;


t_WhiteEnd_TRAIL_green=[];     %TRAIL after white end
t_TRAIL_Death_green=[];        %TRAIL to death
t_TRAIL_Div_F0_green=[];       %after TRAIL
t_TRAIL_WhiteEnd_F1_green=[];  %after TRAIL
t_TRAIL_Div_F1_green=[];
t_TRAIL_WhiteEnd_F2_green=[];
numb=[];

for i=1:length(generation1)
    if generation1(i) == 0 && ~isnan(celldeath1(i))   %cells die in F0
        if whiteendgreenstart1(i) < celldeath1(i) && ~isnan(whiteendgreenstart1(i))   % cells died in SG2M
            t_WhiteEnd_TRAIL_green=[t_WhiteEnd_TRAIL_green, (time(i)-whiteendgreenstart1(i))*0.25];
            t_TRAIL_Death_green=[t_TRAIL_Death_green, (celldeath1(i)-time(i))*0.25];
            t_TRAIL_Div_F0_green=[t_TRAIL_Div_F0_green, nan];
            t_TRAIL_WhiteEnd_F1_green=[t_TRAIL_WhiteEnd_F1_green, nan];
            t_TRAIL_Div_F1_green=[t_TRAIL_Div_F1_green, nan];   
            t_TRAIL_WhiteEnd_F2_green=[t_TRAIL_WhiteEnd_F2_green, nan];
            numb=[numb,cell1(i)];
        else %cells died white (not possible here)
            i;
        end
    elseif generation1(i) == 1 && ~isnan(celldeath1(i))   %cells died in F1
        for m=1:4   
                if round(cell1(i))==round(cell1(i-m)) && generation1(i-m)==0      %i-m was mother
                    if whiteendgreenstart1(i) < celldeath1(i) && ~isnan(whiteendgreenstart1(i))   % cells died in SG2M F1
                        t_WhiteEnd_TRAIL_green=[t_WhiteEnd_TRAIL_green, (time(i)-whiteendgreenstart1(i-m))*0.25];
                        t_TRAIL_Death_green=[t_TRAIL_Death_green, (celldeath1(i)-time(i))*0.25];
                        t_TRAIL_Div_F0_green=[t_TRAIL_Div_F0_green, (celldivision1(i-m)-time(i))*0.25];
                        t_TRAIL_WhiteEnd_F1_green=[t_TRAIL_WhiteEnd_F1_green, (whiteendgreenstart1(i)-time(i))*0.25];
                        t_TRAIL_Div_F1_green=[t_TRAIL_Div_F1_green, nan]; 
                        t_TRAIL_WhiteEnd_F2_green=[t_TRAIL_WhiteEnd_F2_green, nan];
                        numb=[numb,cell1(i)];
                    else    %cells died white F1                        
                        t_WhiteEnd_TRAIL_green=[t_WhiteEnd_TRAIL_green, (time(i)-whiteendgreenstart1(i-m))*0.25];
                        t_TRAIL_Death_green=[t_TRAIL_Death_green, (celldeath1(i)-time(i))*0.25];
                        t_TRAIL_Div_F0_green=[t_TRAIL_Div_F0_green, (celldivision1(i-m)-time(i))*0.25];
                        t_TRAIL_WhiteEnd_F1_green=[t_TRAIL_WhiteEnd_F1_green, nan];
                        t_TRAIL_Div_F1_green=[t_TRAIL_Div_F1_green, nan];
                        t_TRAIL_WhiteEnd_F2_green=[t_TRAIL_WhiteEnd_F2_green, nan];
                        numb=[numb,cell1(i)];
                    end
            end
        end
    elseif generation1(i) == 2 && ~isnan(celldeath1(i))   %cells died in F2
        for m=1:6                                         
            for k=1:6
                    if round(cell1(i))==round(cell1(i-m)) && generation1(i-m)==0 && round(10*cell1(i))==round(10*cell1(i-k)) && generation1(i-k)==1     %i-m was grandmother (f0), i-k was mother (f1)
                        if whiteendgreenstart1(i) < celldeath1(i) && ~isnan(whiteendgreenstart1(i))         % cells died in SG2M F2
                            t_WhiteEnd_TRAIL_green=[t_WhiteEnd_TRAIL_green, (time(i)-whiteendgreenstart1(i-m))*0.25];
                            t_TRAIL_Death_green=[t_TRAIL_Death_green, (celldeath1(i)-time(i))*0.25];
                            t_TRAIL_Div_F0_green=[t_TRAIL_Div_F0_green, (celldivision1(i-m)-time(i))*0.25];
                            t_TRAIL_WhiteEnd_F1_green=[t_TRAIL_WhiteEnd_F1_green, (whiteendgreenstart1(i-k)-time(i))*0.25];
                            t_TRAIL_Div_F1_green=[t_TRAIL_Div_F1_green, (celldivision1(i-k)-time(i))*0.25];   
                            t_TRAIL_WhiteEnd_F2_green=[t_TRAIL_WhiteEnd_F2_green, (whiteendgreenstart1(i)-time(i))*0.25];
                            numb=[numb,cell1(i)];
                        else    %cells died white F1                        
                            t_WhiteEnd_TRAIL_green=[t_WhiteEnd_TRAIL_green, (time(i)-whiteendgreenstart1(i-m))*0.25];
                            t_TRAIL_Death_green=[t_TRAIL_Death_green, (celldeath1(i)-time(i))*0.25];
                            t_TRAIL_Div_F0_green=[t_TRAIL_Div_F0_green, (celldivision1(i-m)-time(i))*0.25];
                            t_TRAIL_WhiteEnd_F1_green=[t_TRAIL_WhiteEnd_F1_green, (whiteendgreenstart1(i-k)-time(i))*0.25];
                            t_TRAIL_Div_F1_green=[t_TRAIL_Div_F1_green, (celldivision1(i-k)-time(i))*0.25];
                            t_TRAIL_WhiteEnd_F2_green=[t_TRAIL_WhiteEnd_F2_green, nan];
                            numb=[numb,cell1(i)];
                        end
                    end
            end
        end
    end    
end



f= [t_WhiteEnd_TRAIL_green;t_TRAIL_Div_F0_green;t_TRAIL_WhiteEnd_F1_green;t_TRAIL_Div_F1_green;t_TRAIL_Death_green;t_TRAIL_WhiteEnd_F2_green];






end

