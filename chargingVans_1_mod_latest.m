%%%% van charging profile for more than one store
%% Charge EVs over night
%% It is assumed that EVs are fully charged at the beginning of each day
%% A slow charger (3 KWh) is used to charger EVs when they are back from their journeys + offset
%% EVs are being charged over night with different chargers
%% Find each store chargers type and number of chargers
%% When there is no more journeys to do, it is assumed that EVs are charging

load branch_idindex_file_test.mat; %% load branch_id indices
branch_idindex;

load branch_id_file.mat branch_id; %% sort branch_id
branch_id;
size_bi = size(branch_id,1);

load sortbranch_id_file.mat sortbranch_id; %% sort branch_id
sortbranch_id;

savedir = 'C:\LocalStore\nt32\Downloads\EvTest1_matlab\output\'; %% new directory

n_vans=[]; %% save branch_id and vansNumber for ev

for id=1:size_bi
    %% inputs
    initial_SoC = 1;
    pack_capacity = 56; %% 75
    energy_per_mile = 0.48;
    charger_rating = 3; %% 12
    journey_safety_margin = 0; %% 0.07 or 0.3*pack_capacity
    HH_energy_transfer = charger_rating/2;
    start_energy = initial_SoC*pack_capacity;
    c_delay = 5; %% time in hour
    cs_offset = c_delay; %% delay between van arrival and charging start (in HH periods)
    van_number = 1; %% fleet size
    max_mpr = 100; %% maximal mileage an EV van do
    
    %% read route information
    %journeys = readtable([ num2str(id) '_newstore_', num2str(sortbranch_id(branch_idindex(id))), '.csv']);
    journeys = readtable([ num2str(id) '_new_store_', num2str(sortbranch_id(branch_idindex(id))), '.csv']);
    %journeys = readtable('1_new_store_101.csv');
    
    T1=table(zeros(size(journeys,1),1), 'VariableNames', {'Energy_Required'}); %% energy required
    T2=table(zeros(size(journeys,1),1), 'VariableNames', {'Energy_Needed'}); %% needed
    journeys = [journeys, T1,T2];
    clear T1;
    clear T2;
    
    numberOfjourneys = size(journeys,1);  %% number of journeys
    dates_j = journeys.Date_Route_Start; %% dates
    dates_j_string = string(dates_j); %% convert dates to string
    formatIn = 'dd/mm/yyyy'; %% date format
    dates_ = datenum(dates_j_string,formatIn); %% format dates into number
    sortDates_ = sort(dates_); %%  sort dates
    datesindex = [1,find(diff(sortDates_))'+1]; %% unique dates index
    datesindex;
    dates = sortDates_(datesindex); %% vector of dates
    numberOfDays = size(datesindex,2); %% number of unique dates
    numberOfRows = 48*numberOfDays; %% each day has 48 HH periods
    
    first_date = min(sortDates_); %% first date on which a journey happens
    last_date = max(sortDates_); %% last date on which a journey happens
    date_size = size(dates,1);  %% total number of days
    
    mpr1 = journeys.Total_MPR*energy_per_mile; %% energy required for each journey column
    mpr = array2table(mpr1); %% converting vector into table
    journeys(:,(end-1)) = mpr; %% energy required for each journey column
    
    eng1 = journeys.Energy_Required*(1+journey_safety_margin); %%  energy needed for each journey column
    eng = array2table(eng1); %% converting vector into table
    journeys(:,end) = eng; %%  energy needed for each journey column
    
    %% HH period per day
    hh = mod([0:(numberOfRows-1)]',48);
    
    %% number of days
    datevector=1:numberOfRows;  %% create dates vector
    for j=0:(numberOfRows-1)
        datevector(j+1)=dates(floor(j/48)+1);
    end
    
    %% ev vans size
    numvans = 200;
    
    %% HH period, dates, time, energy in battery, HH_energy_transfer
    van_energy_profile = zeros(date_size*48,((numvans*2) + 2));
    
    van_energy_profile(:,1) = datevector'; %% dates
    van_energy_profile(:,2) = hh; %% number of HH period
    van_energy_profile(1,3:2:(end-1)) = start_energy*ones(1,numvans); %% adding start_energy to each van Energy column to keep track of the energy
    van_energy_profile(:,4:2:end) = HH_energy_transfer*ones(date_size*48,numvans); %% adding charging rate at all times
    van_energy_profile(48:48:end,4:2:end) = start_energy*ones(date_size,numvans); %% this assume that all vans batteries are fully charge over night
    
    %% This is for fully charging vans on days when here is know journeys
    for d=1:(numberOfDays-1)
        gap=(dates(d+1)-dates(d)-1)*48; %% find days where there is know journeys
        van_energy_profile(48*d,4:2:end) = van_energy_profile(48*d,4:2:end)+gap*HH_energy_transfer*ones(1,numvans); %% charge fully vans
    end
    
    journeys = sortrows(journeys, {'Date_Route_Start','Route_Start'}); %% sort dates and start_time
    
    string_total_time = journeys.Total_Time;  %% Total_Time column
    string_tt = string(string_total_time); %% change total_time into string
    
    string_start_time = journeys.Route_Start;  %% start time column
    string_time_st = string(string_start_time); %% change start time into string
    
    string_end_time = journeys.Route_End; %% end time column
    string_time_end = string(string_end_time); %% change start time into string
    
    matrixOfVans = ones(numberOfRows, numvans);  %% initialising matrixOfVans ('1' = free van, '0' = used van)
    matrixOfVans2 = zeros(numberOfDays, numvans); %% matrix to record journeys Total_MPR
    
    %% minimum number of vans used
    minnumvans = 0; %% ev_vans
    routeids = {[]};  %% initialising the cell entry with an empty vector, ev_vans
    
    dates_i = journeys.Date_Route_Start;
    dates_i_string = string(dates_i);
    formatIn = 'dd/mm/yyyy';
    dates_i_ = datenum(dates_i_string,formatIn); %% format dates into number
    
    rownumST=1;
    
    for i=1:numberOfjourneys  %% loop through all rows of journeys
        
        %% check if journeys are feasible
        if(journeys.Energy_Needed(i) > pack_capacity)
            warning('Insufficient energy, journeys cannot be done:');
            journeys.Energy_Needed(i);
            i;
        end
        
        dateSTindex = find(dates==dates_i_(i))-1; %% get nonzeroes elements of the matrix
        
        %% change time into HH: HH = 2*(hour + (minute/60))
        total_time = ceil(2*(hour(string_tt(i)) + (minute(string_tt(i))/60))); %% time into number(minute
        
        rownumSTold=rownumST;
        timeST = floor(2*(hour(string_time_st(i)) + (minute(string_time_st(i))/60))); %% time into number(minute)
        rownumST = dateSTindex*48+timeST+1 ;  %% A day start at 00 and ends at 47
        
        %% update energy matrix from rownumSTold to rownumST
        if rownumST > rownumSTold
            
            for rownums=rownumSTold:(rownumST-1)
                van_energy_profile(rownums+1,3:2:(end-1)) = min(pack_capacity, van_energy_profile(rownums,3:2:(end-1))+van_energy_profile(rownums,4:2:end)); %% energy before journeys start
            end
        end
        
        timeEND = floor(2*(hour(string_time_end(i)) + (minute(string_time_end(i))/60))); %% time into number(minute)
        rownumEND = dateSTindex*48+timeEND+1;  %% A day start at 00 and ends at 47
        
        %% Allocate journeys
        if (journeys.Total_MPR(i) < max_mpr) %% ev_vans
            vantempnum = find(matrixOfVans(rownumST,:)); %% index of free vans
            sizetemp = length(vantempnum);
            
            %% check if there are enough vans to do the journeys
            if(sizetemp==0)
                warning("not enough vans ...!!! ");
                break
            end
            
            mytest=0;
            idx=0; %% idx = 0 pick the first free van
            
            while(mytest==0)
                idx=idx+1; %% pick the next free van
                vannum = vantempnum(idx); %% van number
                
                %% Check if there is enough energy to do next journey and if total_MPR < max_mpr
                if( (van_energy_profile(rownumST, 2*vannum+1)>journeys.Total_MPR(i)*energy_per_mile) & (matrixOfVans2(dateSTindex+1, vannum)+journeys.Total_MPR(i) < max_mpr) )
                    mytest = 1;
                end
            end
            
            matrixOfVans(rownumST:rownumEND,vannum)=zeros(rownumEND-rownumST+1,1); %% updating matrixOfVans
            matrixOfVans2(dateSTindex+1,vannum)=matrixOfVans2(dateSTindex+1,vannum)+journeys.Total_MPR(i);
            van_energy_profile(rownumST:rownumEND,2+2*vannum)=-journeys.Total_MPR(i)*energy_per_mile/(-rownumST+rownumEND+1)*ones(rownumEND-rownumST+1,1);
        end
    end
    
    while(nnz(matrixOfVans(:,minnumvans+1))<numberOfRows && minnumvans<numvans-1)
            minnumvans=minnumvans+1;
    end
        
        %% print the minimum number of vans used
        vansNumber = minnumvans;
        vansNumber;

 %% this is to check if all the used vans have a charger assigned to them
    n_vans=[n_vans;branch_id(id), vansNumber]; %% append
    file_name = fopen('used_vans.csv', 'wt'); %% save 'a' in csv file
    fprintf(file_name, '%s,%s\n', 'Branch_id', 'ev_vansNumber');
    fprintf(file_name, '%g,%g\n', n_vans');
    fclose(file_name);
    
    %% Copy file from Matlab work space to a given folder
    copyfile('used_vans.csv', savedir);

    
    %% update the van energy profile matrix
    for rownums=rownumST:(48*date_size-1)
        van_energy_profile(rownums+1,3:2:(end-1)) = min(pack_capacity, van_energy_profile(rownums,3:2:(end-1))+van_energy_profile(rownums,4:2:end));
    end
    
    %% updating pack_capacity and power flow columns
    for rownums=1:(48*date_size-1)
        van_energy_profile(rownums,4:2:end) = van_energy_profile(rownums+1,3:2:(end-1))-van_energy_profile(rownums,3:2:(end-1));
    end
    
    %% put in the last row of power flow column how much energy is needed for the van to be fully charged
    van_energy_profile(48*date_size,4:2:end) = pack_capacity-van_energy_profile(48*date_size,3:2:(end-1)); %% energy before journeys start
    
    T = array2table(van_energy_profile);
    writetable(T, [num2str(id) '_this_output_' num2str(sortbranch_id(branch_idindex(id))) '.csv'], 'Delimiter',','); %%
    writetable(T, [savedir num2str(id) '_this_output_' num2str(sortbranch_id(branch_idindex(id))) '.csv'], 'Delimiter',','); %% save output in a different folder
    
    %% This is for charging over night and finding chargers type
    charger_required = zeros(date_size, vansNumber); %% matrix to records chargers type
    
    for currentvan=1:vansNumber %% for loop for number of vans used
        returnvector=1+find(diff(matrixOfVans(1:48, currentvan))==1); %% time at which a van is return from a journey
        
        for k=2:(date_size)
            k;
            startvector =1+find(diff(matrixOfVans(((k-1)*48+1):(k*48), currentvan))==-1); %% time at which a van start it journey form day 2
            
            if min(size(startvector))>0 & min(size(returnvector))>0
                numtimeslots = (48*(dates(k)-dates(k-1))-1)+startvector(1)-returnvector(end); %% number of charging slots available to charge over night
                1+2*currentvan; %% energy in the pack column for each van
                van_energy_profile(returnvector(end), 1+2*currentvan); %% energy in the battery before charging start
                
                %% determine charger type
                charger_required(k-1, currentvan) = (pack_capacity - van_energy_profile(returnvector(end), 1+2*currentvan)) / numtimeslots; % required charging rate per 30 minutes (HH)
            else
                if min(size(returnvector))>0
                    numtimeslots = (48*((dates(k)-dates(k-1))+1)-1)-returnvector(end); %% number of HH availble for charging
                    
                    %% determine chargertype
                    charger_required(k-1, currentvan) = (pack_capacity - van_energy_profile(returnvector(end), 1+2*currentvan)) / numtimeslots; % required charging rate per 30 minutes (HH)
                end
            end
            returnvector=1+find(diff(matrixOfVans(((k-1)*48+1):(k*48), currentvan))==1); %% time at which vans are back after journeys
        end
    end
    
    %% record the number of charger required for each van used
    charger_rates = [3,7,11,22]/2; % available rates per 30 minutes (HH)
    [ordered_charger_rates, permute] = sort(charger_rates); % order according to rate. permute for order(1,2,3,4)
    charger_numtypes = length(charger_rates); %% number of different chargers (or chargers type)
    charger_collection = zeros(1,vansNumber); %% row vector to find chargers to be used
    maxcharge = max(sort(charger_required,2)); %% sort 'charger_required' by row then find the max
    
    for k=1:vansNumber
        i=1;
        maxcharge(k);
        while ((maxcharge(k) > ordered_charger_rates(i)) & (i< charger_numtypes))
            i=i+1;
        end
        
        if maxcharge(k) > ordered_charger_rates(i)
            % display error
            warning(" No charger ...!!! ");
            break
        end
        
        charger_collection(k) = ordered_charger_rates(i); %% row matrix of chargers used by each van
    end
    
    charger_collection; %% list of chargers used
    total_chargers = size(charger_collection,2);
    
    number_of_charger_3 = sum(charger_collection(1,:)==charger_rates(1));
    number_of_charger_7 = sum(charger_collection(1,:)==charger_rates(2));
    number_of_charger_11 = sum(charger_collection(1,:)==charger_rates(3));
    number_of_charger_22 = sum(charger_collection(1,:)==charger_rates(4));
    
    %% charger_3, charger_7, charger_11, charger_22, total number of chargers
    chargers_num = [number_of_charger_3,number_of_charger_7,number_of_charger_11,number_of_charger_22,total_chargers];
    n_chargers = array2table(chargers_num); %% change matrix into table
    col_names = {'charger_3','charger_7','charger_11','charger_22','total_chargers'}; %% table headers
    n_chargers.Properties.VariableNames = col_names; %% table headers
    chargers = n_chargers
    
    writetable(chargers, [num2str(id) 'chargers_type_' num2str(sortbranch_id(branch_idindex(id))) '.csv'], 'Delimiter',','); %%
    writetable(chargers, [savedir num2str(id) '_chargers_type_' num2str(sortbranch_id(branch_idindex(id))) '.csv'], 'Delimiter',','); %% save output in a different folder
    
end




