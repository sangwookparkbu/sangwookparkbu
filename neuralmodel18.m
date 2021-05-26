clear
clc
close all
% all we have to do at the end 1. 
%%dx/dt=-ax+I type equations modeled 
dT=0.005;
timedifference = 10; % % different processing time for intraconnection and interconnection
DT = dT*timedifference;
% speed = 8; %scalar value of speed 1 neuron unit/ 1 s
neuronarray = 200; % neuron number x y has to be 200 for zollner illusion
time = 30;
% time = floor(neuronarray/speed);
% dtime=time/dT;
% k = dtime/time;
% dspeed=speed/k; % 1 neuron unit/5 ms
Ie = zeros(time,neuronarray*2,neuronarray*2);
Ie=shiftdim(Ie,1);
Iix(:,:,:) = Ie(:,:,:).*0;
localtheta = 20;
orientationtheta = 90; % keep it 90
c = 1:30:neuronarray;
% aa=a*pi./180;
f = 0.5; % frequency of local contextual lines. Less magnitude of f means more frequency. Also, lower f means shorter local line
% relationship between frequency of local line and local line length: local line length = frequency*size(localine,1); 
% too short of frequency will not work for high angle

% Ie(neuronarray,neuronarray,1)=1; % for checking receptive field

 % Ie(:,neuronarray,:)=1; % for checking
 % %  Ie(:,:,:)=imrotate(Ie(:,:,:),-20,'bilinear','crop');
% Ie(neuronarray,:,:)=1;
% sortbox = sum(Ie,3);
% [sortx, sorty]=ind2sub(size(sortbox),find(sortbox>0.8.*max(sortbox,[],1:2))); % % this is to look at the individual hypercolumn responses later
% [sortedx,sorted]=sort(sortx);
% sortedy = sorty(sorted);
% vectorbox = sum(Ie,3);
% [vecx, vecy]=ind2sub(size(vectorbox),find(vectorbox>0.6.*max(vectorbox,[],1:2))); % % this is to sort out relevant vector plotting later
% [vectorx,sortedyy]=sort(vecx);
% vectory = vecy(sortedyy);

% Ie(neuronarray:neuronarray+10,neuronarray,:)=1; % for testing corner detection
% Ie(neuronarray,neuronarray:neuronarray+10,:)=1; %

% % % % ponzo illusion
% localine = zeros(150,150,size(Ie,3)); % for ponzo illusion horizontal
% localine(size(localine,1)*0.5,:,:) = 1;
% localinee = localine;
% localine(:,:,:) = imrotate(localine,-localtheta,'bilinear','crop');
% localinee(:,:,:) = imrotate(localinee,localtheta,'bilinear','crop');
% Ie(neuronarray*1+1:neuronarray*1.75,neuronarray*0.625+1:neuronarray*1.375,:) = localine(:,:,:);
% Ie(neuronarray*0.25+1:neuronarray*1,neuronarray*0.625+1:neuronarray*1.375,:)=localinee(:,:,:);
% Ie(326:end,:,:)=[];
%  Ie(1:neuronarray*0.375,:,:)=[];
%  Ie(:,301:end,:)=[];
%  Ie(:,1:neuronarray*0.5,:)=[];
% Ie(101:150,:,:)=[];
% Ie(81:120,169,:)=1;
% Ie(81:120,30,:)=1;


%  %%ponzo illusion
% localine = zeros(150,150,size(Ie,3)); % for ponzo illusion
% localine(:,size(localine,2)*0.5,:) = 1;
% localinee = localine;
% localine(:,:,:) = imrotate(localine,-localtheta,'bilinear','crop');
% localinee(:,:,:) = imrotate(localinee,localtheta,'bilinear','crop');
% Ie(neuronarray*0.625+1:neuronarray*1.375,neuronarray*1+1:neuronarray*1.75,:) = localine(:,:,:);
% Ie(neuronarray*0.625+1:neuronarray*1.375,neuronarray*0.25+1:neuronarray*1,:)=localinee(:,:,:);
% Ie(:,326:end,:)=[];
%  Ie(:,1:neuronarray*0.375,:)=[];
%  Ie(301:end,:,:)=[];
%  Ie(1:neuronarray*0.5,:,:)=[];
% Ie(:,101:150,:)=[];
% Ie(169,81:120,:)=1;
% Ie(35,81:120,:)=1;
% 
% % %poggendorff illusion
Ip = Ie;
Ip2 = Ie;
localine = zeros(200,200,size(Ie,3));
localine(:,size(localine,2)*0.5,:) = 1;
 localine(:,:,:) = imrotate(localine,localtheta,'bilinear','crop');
 localinee = localine;
 localinee(101:end,:,:)=0;
Ip(neuronarray*0.5+1:neuronarray*1.5,neuronarray*0.5+1:neuronarray*1.5,:)=localine(:,:,:);
% Ip2(neuronarray*0.5+1:neuronarray*1.5,neuronarray*0.5+1-4:neuronarray*1.5-4,:)=localinee(:,:,:);
Ie = Ip + Ip2;
Ie(:,326:end,:)=[];
Ie(:,1:75,:)=[];
Ie(326:end,:,:)=[];
Ie(1:75,:,:)=[];
Ie(:,116:134,:)=0;
sortbox = sum(Ie,3);
[sortx, sorty]=ind2sub(size(sortbox),find(sortbox>0.8.*max(sortbox,[],1:2))); % % this is to look at the individual hypercolumn responses later
[sortedx,sorted]=sort(sortx);
sortedy = sorty(sorted);
Ie(1:end,115,:)=1;
Ie(1:end,135,:)=1;
vectorbox = sum(Ie,3);
[vecx, vecy]=ind2sub(size(vectorbox),find(vectorbox>0.6.*max(vectorbox,[],1:2))); % % this is to sort out relevant vector plotting later
[vectorx,sortedyy]=sort(vecx);
vectory = vecy(sortedyy);
% Ie(:,116:134,:)=0;

% % % % this is the atypical poggendorff illusion for demonstration purpose
% Ip = Ie;
% Ip2 = Ie;
% localine = zeros(200,200,size(Ie,3));
% localine(:,size(localine,2)*0.5,:) = 1;
%  localine(:,:,:) = imrotate(localine,localtheta,'bilinear','crop');
%  localinee = localine;
%  localinee(101:end,:,:)=0;
% Ip(neuronarray*0.5+1:neuronarray*1.5,neuronarray*0.5+1:neuronarray*1.5,:)=localine(:,:,:);
% % Ip2(neuronarray*0.5+1:neuronarray*1.5,neuronarray*0.5+1-4:neuronarray*1.5-4,:)=localinee(:,:,:);
% Ie = Ip + Ip2;
% Ie(:,326:end,:)=[];
% Ie(:,1:75,:)=[];
% Ie(326:end,:,:)=[];
% Ie(1:75,:,:)=[];
% sortbox = sum(Ie,3);
% [sortx, sorty]=ind2sub(size(sortbox),find(sortbox>0)); % % this is to look at the individual hypercolumn responses later
% [sortedx,sorted]=sort(sortx);
% sortedy = sorty(sorted);
% Ie(1:end,115,:)=1;
% Ie(1:end,135,:)=1;



% % % % % % % % rod and frame
% localbox = zeros(150,150,size(Ie,3));
% localbox(1:end,1,:)=1;
% localbox(1:end,end,:)=1;
% localbox(1,1:end,:)=1;
% localbox(end,1:end,:)=1;
% ab = (size(Ie,1)-size(localbox,1))*0.5;
% Ie(ab+1:ab+size(localbox,1),ab+1:ab+size(localbox,1),:)=localbox(:,:,:);
% %Ie(:,:,:) = imrotate(Ie,-localtheta,'bilinear','crop');
% Ie(301:end,:,:)=[];
% Ie(1:100,:,:)=[];
% Ie(:,301:end,:)=[];
% Ie(:,1:100,:)=[];
% Ie(neuronarray*0.5-55+1:neuronarray*0.5+55,neuronarray*0.5,:)=1;

% localine = zeros(80,80,size(Ie,3)); % for zollner illusion horizontal
% localine(size(localine,1)*0.5,:,:) = 1;
% localinee = localine;
% localine(:,:,:) = imrotate(localine,-localtheta,'bilinear','crop');
% localinee(:,:,:) = imrotate(localinee,localtheta,'bilinear','crop');
% % c = 1:neuronarray*2/size(localine,1);
% for a = 1:f:15 % frequency of the local contextual line
%     Ie(neuronarray*1.3+1:neuronarray*1.7,a*size(localine,1)-size(localine,1)+1:a*size(localine,1),:) = localine(:,:,:);
%     Ie(neuronarray*0.9+1:neuronarray*1.3,size(localinee,1)*a-size(localinee,1)+1:size(localinee,1)*a,:)=localinee(:,:,:);
%     Ie(neuronarray*0.5+1:neuronarray*0.9,size(localine,1)*a-size(localine,1)+1:size(localine,1)*a,:)=localine(:,:,:);
% end
% % Ie(30,:,:)=1;
% Ie(:,346:end,:)=[];
% Ie(:,1:neuronarray*0.475,:)=[];
% Ie(346:end,:,:)=[];
% Ie(1:neuronarray*0.475,:,:)=[];
% Ie(round(((neuronarray*1.025-40*sind(localtheta)))+(f*80*sind(localtheta)*0.5)),:,:)=1;
% Ie(round(((neuronarray*1.025-40*sind(-localtheta)))+(f*80*sind(-localtheta)*0.5))-80,:,:)=1;
% Ie(round(((neuronarray*1.025-40*sind(localtheta)))+(f*80*sind(localtheta)*0.5))-160,:,:)=1;


% localine = zeros(80,80,size(Ie,3)); % for zollner illusion vertical
% localine(size(localine,1)*0.5,:,:) = 1;
% localinee = localine;
% localine(:,:,:) = imrotate(localine,90+localtheta,'bilinear','crop');
% localinee(:,:,:) = imrotate(localinee,90-localtheta,'bilinear','crop');
% % c = 1:neuronarray*2/size(localine,1);
% for a = 1:f:15 % frequency of the local contextual line
%     Ie(a*size(localine,1)-size(localine,1)+1:a*size(localine,1),neuronarray*1.3+1:neuronarray*1.7,:) = localine(:,:,:);
%     Ie(size(localinee,1)*a-size(localinee,1)+1:size(localinee,1)*a,neuronarray*0.9+1:neuronarray*1.3,:)=localinee(:,:,:);
%     Ie(size(localine,1)*a-size(localine,1)+1:size(localine,1)*a,neuronarray*0.5+1:neuronarray*0.9,:)=localine(:,:,:);
% end
% % Ie(30,:,:)=1;
% Ie(346:end,:,:)=[];
% Ie(1:neuronarray*0.475,:,:)=[];
% Ie(:,346:end,:)=[];
% Ie(:,1:neuronarray*0.475,:)=[];
% Iix(346:end,:,:)=[];
% Iix(1:neuronarray*0.475,:,:)=[];
% Iix(:,346:end,:)=[];
% Iix(:,1:neuronarray*0.475,:)=[];
% Iix(:,round(((neuronarray*1.025-40*sind(localtheta)))+(f*80*sind(localtheta)*0.5)),:)=1;
% Iix(:,round(((neuronarray*1.025-40*sind(-localtheta)))+(f*80*sind(-localtheta)*0.5))-80,:)=1;
% Iix(:,round(((neuronarray*1.025-40*sind(localtheta)))+(f*80*sind(localtheta)*0.5))-160,:)=1;
% sortbox = sum(Iix,3);
% [sortx, sorty]=ind2sub(size(sortbox),find(sortbox>0)); % % this is to look at the individual hypercolumn responses later
% [sortedx,sorted]=sort(sortx);
% sortedy = sorty(sorted);
% Ie = Ie + Iix;
% Ie(:,50:78,:)=[];
% vectorbox = sum(Ie,3);
% [vecx, vecy]=ind2sub(size(vectorbox),find(vectorbox>0.3.*max(vectorbox,[],1:2))); % % this is to look at the individual hypercolumn responses later
% [vectorx,sortedyy]=sort(vecx);
% vectory = vecy(sortedyy);



% for a = 1:size(Ie,1)
% for b = 1:size(Ie,2)
%     if Ie(a,b,:)<-0.5
%         Ie(a,b,:)=-1;
%     end
% end
% end

% imagesc(Ie(:,:,10));
% colormap(flipud(gray(10)));

% localine = zeros(80,80,size(Ie,3)); % for zollner illusion vertical one line zollner vertical
% localine(size(localine,1)*0.5,:,:) = 1;
% localinee = localine;
% localine(:,:,:) = imrotate(localine,90+localtheta,'bilinear','crop');
% localinee(:,:,:) = imrotate(localinee,90-localtheta,'bilinear','crop');
% % c = 1:neuronarray*2/size(localine,1);
% for a = 1:f:15 % frequency of the local contextual line
%   %  Ie(a*size(localine,1)-size(localine,1)+1:a*size(localine,1),neuronarray*1.3+1:neuronarray*1.7,:) = localine(:,:,:);
%     Ie(size(localinee,1)*a-size(localinee,1)+1:size(localinee,1)*a,neuronarray*0.9+1:neuronarray*1.3,:)=localinee(:,:,:);
%   %  Ie(size(localine,1)*a-size(localine,1)+1:size(localine,1)*a,neuronarray*0.5+1:neuronarray*0.9,:)=localine(:,:,:);
% end
% % Ie(30,:,:)=1;
% Ie(346:end,:,:)=[];
% Ie(1:neuronarray*0.475,:,:)=[];
% Ie(:,346:end,:)=[];
% Ie(:,1:neuronarray*0.475,:)=[];
% Iix(346:end,:,:)=[];
% Iix(1:neuronarray*0.475,:,:)=[];
% Iix(:,346:end,:)=[];
% Iix(:,1:neuronarray*0.475,:)=[];
% %Ie(:,round(((neuronarray*1.025-40*sind(localtheta)))+(f*80*sind(localtheta)*0.5)),:)=1;
% Iix(:,round(((neuronarray*1.025-40*sind(-localtheta)))+(f*80*sind(-localtheta)*0.5))-80,:)=1;
% sidesortbox = sum(Ie,3);
% [sidesortx, sidesorty]=ind2sub(size(sidesortbox),find(sidesortbox>0)); % % this is to look at the individual hypercolumn responses later
% [sidesortedx,sidesorted]=sort(sidesortx);
% sidesortedy = sidesorty(sidesorted);
% sortbox = sum(Iix,3);
% [sortx, sorty]=ind2sub(size(sortbox),find(sortbox>0)); % % this is to look at the individual hypercolumn responses later
% [sortedx,sorted]=sort(sortx);
% sortedy = sorty(sorted);
% Ie(:,:,:) = Ie + Iix;
% vectorbox = sum(Ie,3);
% [vecx, vecy]=ind2sub(size(vectorbox),find(vectorbox>0.3.*max(vectorbox,[],1:2))); % % this is to look at the individual hypercolumn responses later
% [vectorx,sortedyy]=sort(vecx);
% vectory = vecy(sortedyy);


%Ie(:,round(((neuronarray*1.025-40*sind(localtheta)))+(f*80*sind(localtheta)*0.5))-160,:)=1;





% Ie(:,size(Ie,2)*0.5,:) = 1; % % % this is for investigating sub additivity of normalization
% Ie(:,:,:) = imrotate(Ie,localtheta,'bilinear','crop');
% % Ie(30,:,:)=1;
% Ie(:,301:end,:)=[];
% Ie(:,1:neuronarray*0.5,:)=[];
% Ie(301:end,:,:)=[];
% Ie(1:neuronarray*0.5,:,:)=[];
% Ie(:,size(Ie,2)*0.5,:)=1;


% for a = 1:neuronarray*2-1 % for 2 crossing lines this will be used for 
% %Ie(20+a,1+c+a,:)=1;
% 
% 
%  Iix(neuronarray,end-a,:)=1;
% end
% Iix(:,:,:) = imrotate(Iix,localtheta,'bilinear','crop');
% % Ie(30,:,:)=1;
% Iix(:,301:end,:)=[];
% Iix(:,1:neuronarray*0.5,:)=[];
% Iix(301:end,:,:)=[];
% Iix(1:neuronarray*0.5,:,:)=[];
% Iix(neuronarray*0.5,:,:)=1;

%   for a = 1:30 % for 2 line zollner illusion 
%        Ie(150+a,1+c+a,:)=1;
% 
% 
%       Ie(250-a,1+c+a,:)=1;
%   end
% Ie(:,201:end,:)=[];
% % Ie(:,1:neuronarray*0.5,:)=[];
% Ie(301:end,:,:)=[];
% Ie(1:neuronarray*0.5,:,:)=[];
% Ie(65,:,:)=1;
% Ie(135,:,:)=1;




 %  n2=size(Ie,2)*size(Ie,3); % number of neurons
%     dT=0.005;   %differntial dT set at 0.05 seconds, can be set at any value of diff time
 m=max(size(Ie,3))*dT;    %gives max T a connection to the amount of inputs making sure it is the proper number to make #t=#I 
                            %works as long as all inputs have the same
                            %dimensions, if not different neurons can have
                            %different values
T=m;    % max 'T' set at m seconds
tt=0:dT:T;    %setting values for t
n=max(tt/dT);       %connects t and i

 
%need to add weights


RGCsige = 1;
LGNsige = 5;
V1sige = 1;
vhratio = 9; % V1 RF vertical length:horizontal length ratio. increasing this would increase orientation sensitivity

eiratio = 2.1; % the ratio between excitatory field and inhibitory field in v1. evidence suggests that low level processing is intact in asd
eiratiov2 = 2.1; % inhibitory/excitatory ratio for v2
iceiratio = 2.1; % inhibitory/excitatory ratio for interaction field
sizeconstant = 1; % increasing this would increase the receptive field size of V2
RGCe = floor(RGCsige*eiratio*7.4);
LGNe = floor(LGNsige*eiratio*7.4);
V1e = floor(vhratio*7.4*V1sige);
V2e = floor(vhratio*7.4*(V1sige+sizeconstant));
xsigmae = 3;
ysigmae = 3;
zsigmae = 3;

if mod(round(LGNe),2) == 0 % % this is important for local interaction balance
    j = 1;
elseif mod(round(LGNe),2) == 1
    j = 0;
end

if mod(round(RGCe),2) == 0 % % this is important for local interaction balance
    r = 1;
elseif mod(round(RGCe),2) == 1
    r = 0;
end

if mod(round(V1e),2) == 0 % % this is important for local interaction balance
    m = 1;
elseif mod(round(V1e),2) == 1
    m = 0;
end

if mod(round(V2e),2) == 0 % % this is important for local interaction balance
    w = 1;
elseif mod(round(V2e),2) == 1
    w = 0;
end

RGCee=(fspecial('gaussian',[RGCe+r RGCe+r], RGCsige)); % excitation normal gaussian 
RGCii=(fspecial('gaussian',[RGCe+r RGCe+r], RGCsige*eiratio)); % inhibition 

% LGNee=(fspecial('gaussian',[LGNe+j LGNe+j], LGNsige)); % excitation normal gaussian 
% LGNii=(fspecial('gaussian',[LGNe+j LGNe+j], LGNsige*eiratio)); % inhibition 

LGNee=(fspecial('gaussian',[121 121], LGNsige)); % excitation normal gaussian 
LGNii=(fspecial('gaussian',[121 121], LGNsige*eiratio)); % inhibition 

%  WV1iv=fspecial('gaussian',[V1e V1e], (4)); % concentric receptive field
%  WV1ih=fspecial('gaussian',[V1e V1e], (4)*eiratio); % 


WV1ev=fspecial('gaussian',[V1e+m V1e+m], V1sige);  
WV1eh=fspecial('gaussian',[V1e+m V1e+m], (V1sige)*vhratio); 

WV1iv=fspecial('gaussian',[V1e+m V1e+m], (V1sige)*eiratio); % 
WV1ih=fspecial('gaussian',[V1e+m V1e+m], (V1sige)*vhratio); % 

%WV1ih=fspecial('gaussian',[V1e V1e], [1 9]);


WV2ev=fspecial('gaussian',[V2e+w V2e+w], V1sige+sizeconstant);  
WV2eh=fspecial('gaussian',[V2e+w V2e+w], (V1sige+sizeconstant)*vhratio); 

WV2iv=fspecial('gaussian',[V2e+w V2e+w], (V1sige+sizeconstant)*(eiratiov2)); % V2 receptive field is more elliptical than V1
WV2ih=fspecial('gaussian',[V2e+w V2e+w], (V1sige+sizeconstant)*vhratio);

if mod(round(zsigmae*6),2) == 0 % % this is important for local interaction balance
    a = 1;
elseif mod(round(zsigmae*6),2) == 1
    a = 0;
end

if mod(round(xsigmae*6),2) == 0 % % this is important for local interaction balance
    c= 1;
elseif mod(round(xsigmae*6),2) == 1
    c = 0;
end

if mod(round(zsigmae*iceiratio*6),2) == 0
    b = 1;
elseif mod(round(zsigmae*iceiratio*6),2) == 1
    b = 0;
end

if mod(round(xsigmae*iceiratio*6),2) == 0
    d = 1;
elseif mod(round(xsigmae*iceiratio*6),2) == 1
    d = 0;
end

Interactione = fspecial3('gaussian',[round(xsigmae*6)+c,round(ysigmae*6)+c,round(zsigmae*6)+a],[xsigmae ysigmae zsigmae]);
Interactioni = fspecial3('gaussian',[round(xsigmae*iceiratio*6)+d, round(ysigmae*iceiratio*6)+d,round(zsigmae*iceiratio*6+b)],[xsigmae*iceiratio ysigmae*iceiratio  zsigmae*iceiratio]);

%  WV1ev=fspecial('gaussian',[V1e V1e], 1+sizeconstant);  
% WV1eh=fspecial('gaussian',[V1e V1e], (1+sizeconstant)*vhratio); 
%  
%  WV1iv=fspecial('gaussian',[V1e V1e], (1+sizeconstant)*eiratio); % for eggshaped receptive field
%  WV1ih=fspecial('gaussian',[V1e V1e], (1+sizeconstant)*vhratio*eiratio); % 

V1horizontale = WV1ev*WV1eh;
V1horizontali = WV1iv*WV1ih;

V1horizontalee = WV1ev*WV1eh;
V1horizontalii = WV1iv*WV1ih;

V2horizontale = WV2ev*WV2eh;
V2horizontali = WV2iv*WV2ih;

V2horizontalee = WV2ev*WV2eh;
V2horizontalii = WV2iv*WV2ih;

V2verticale = imrotate(V2horizontalee,orientationtheta,'bilinear','crop');
V2verticali = imrotate(V2horizontalii,orientationtheta,'bilinear','crop');

V210i = imrotate(V2horizontalii,-10,'bilinear','crop');
V210e = imrotate(V2horizontalee,-10,'bilinear','crop');
V220i = imrotate(V2horizontalii,-20,'bilinear','crop');
V220e = imrotate(V2horizontalee,-20,'bilinear','crop');

V230e = imrotate(V2horizontalee,-30,'bilinear','crop');
V230i = imrotate(V2horizontalii,-30,'bilinear','crop');

V240e = imrotate(V2horizontalee,-40,'bilinear','crop');
V240i = imrotate(V2horizontalii,-40,'bilinear','crop');
V250e = imrotate(V2horizontalee,-50,'bilinear','crop');
V250i = imrotate(V2horizontalii,-50,'bilinear','crop');

V260e = imrotate(V2horizontalee,-60,'bilinear','crop');
V260i = imrotate(V2horizontalii,-60,'bilinear','crop');

V270e = imrotate(V2horizontalee,-70,'bilinear','crop');
V270i = imrotate(V2horizontalii,-70,'bilinear','crop');

V280e = imrotate(V2horizontalee,-80,'bilinear','crop');
V280i = imrotate(V2horizontalii,-80,'bilinear','crop');

V2100e = imrotate(V2horizontalee,-100,'bilinear','crop');
V2100i = imrotate(V2horizontalii,-100,'bilinear','crop');

V2110e = imrotate(V2horizontalee,-110,'bilinear','crop');
V2110i = imrotate(V2horizontalii,-110,'bilinear','crop');

V2120e = imrotate(V2horizontalee,-120,'bilinear','crop');
V2120i = imrotate(V2horizontalii,-120,'bilinear','crop');

V2130e = imrotate(V2horizontalee,-130,'bilinear','crop');
V2130i = imrotate(V2horizontalii,-130,'bilinear','crop');

V2140e = imrotate(V2horizontalee,-140,'bilinear','crop');
V2140i = imrotate(V2horizontalii,-140,'bilinear','crop');

V2150e = imrotate(V2horizontalee,-150,'bilinear','crop');
V2150i = imrotate(V2horizontalii,-150,'bilinear','crop');

V2160e = imrotate(V2horizontalee,-160,'bilinear','crop');
V2160i = imrotate(V2horizontalii,-160,'bilinear','crop');

V2170e = imrotate(V2horizontalee,-170,'bilinear','crop');
V2170i = imrotate(V2horizontalii,-170,'bilinear','crop');


V2verticalee = imrotate(V2horizontalee,orientationtheta,'bilinear','crop');
V2verticalii = imrotate(V2horizontalii,orientationtheta,'bilinear','crop');

V230ee = imrotate(V2horizontalee,-30,'bilinear','crop');
V230ii = imrotate(V2horizontalii,-30,'bilinear','crop');

V260ee = imrotate(V2horizontalee,-60,'bilinear','crop');
V260ii = imrotate(V2horizontalii,-60,'bilinear','crop');


V2120ee = imrotate(V2horizontalee,-120,'bilinear','crop');
V2120ii = imrotate(V2horizontalii,-120,'bilinear','crop');
V2135ee = imrotate(V2horizontalee,-135,'bilinear','crop');
V2135ii = imrotate(V2horizontalii,-135,'bilinear','crop');
V2150ee = imrotate(V2horizontalee,-150,'bilinear','crop');
V2150ii = imrotate(V2horizontalii,-150,'bilinear','crop');
V2165ee = imrotate(V2horizontalee,-165,'bilinear','crop');
V2165ii = imrotate(V2horizontalii,-165,'bilinear','crop');

V1verticale = imrotate(V1horizontalee,orientationtheta,'bilinear','crop');
V1verticali = imrotate(V1horizontalii,orientationtheta,'bilinear','crop');

V110i = imrotate(V1horizontalii,-10,'bilinear','crop');
V110e = imrotate(V1horizontalee,-10,'bilinear','crop');
V120i = imrotate(V1horizontalii,-20,'bilinear','crop');
V120e = imrotate(V1horizontalee,-20,'bilinear','crop');

V130e = imrotate(V1horizontalee,-30,'bilinear','crop');
V130i = imrotate(V1horizontalii,-30,'bilinear','crop');

V140e = imrotate(V1horizontalee,-40,'bilinear','crop');
V140i = imrotate(V1horizontalii,-40,'bilinear','crop');

V150e = imrotate(V1horizontalee,-50,'bilinear','crop');
V150i = imrotate(V1horizontalii,-50,'bilinear','crop');

V160e = imrotate(V1horizontalee,-60,'bilinear','crop');
V160i = imrotate(V1horizontalii,-60,'bilinear','crop');

V170e = imrotate(V1horizontalee,-70,'bilinear','crop');
V170i = imrotate(V1horizontalii,-70,'bilinear','crop');

V180e = imrotate(V1horizontalee,-80,'bilinear','crop');
V180i = imrotate(V1horizontalii,-80,'bilinear','crop');


V1100e = imrotate(V1horizontalee,-100,'bilinear','crop');
V1100i = imrotate(V1horizontalii,-100,'bilinear','crop');

V1110e = imrotate(V1horizontalee,-110,'bilinear','crop');
V1110i = imrotate(V1horizontalii,-110,'bilinear','crop');

V1120e = imrotate(V1horizontalee,-120,'bilinear','crop');
V1120i = imrotate(V1horizontalii,-120,'bilinear','crop');
V1130e = imrotate(V1horizontalee,-130,'bilinear','crop');
V1130i = imrotate(V1horizontalii,-130,'bilinear','crop');

V1140e = imrotate(V1horizontalee,-140,'bilinear','crop');
V1140i = imrotate(V1horizontalii,-140,'bilinear','crop');

V1150e = imrotate(V1horizontalee,-150,'bilinear','crop');
V1150i = imrotate(V1horizontalii,-150,'bilinear','crop');

V1160e = imrotate(V1horizontalee,-160,'bilinear','crop');
V1160i = imrotate(V1horizontalii,-160,'bilinear','crop');

V1170e = imrotate(V1horizontalee,-170,'bilinear','crop');
V1170i = imrotate(V1horizontalii,-170,'bilinear','crop');




V1verticalee = imrotate(V1horizontalee,orientationtheta,'bilinear','crop');
V1verticalii = imrotate(V1horizontalii,orientationtheta,'bilinear','crop');
V115ee = imrotate(V1horizontalee,-15,'bilinear','crop');
V115ii = imrotate(V1horizontalii,-15,'bilinear','crop');
V130ee = imrotate(V1horizontalee,-30,'bilinear','crop');
V130ii = imrotate(V1horizontalii,-30,'bilinear','crop');
V145ee = imrotate(V1horizontalee,-45,'bilinear','crop');
V145ii = imrotate(V1horizontalii,-45,'bilinear','crop');
V160ee = imrotate(V1horizontalee,-60,'bilinear','crop');
V160ii = imrotate(V1horizontalii,-60,'bilinear','crop');
V175ee = imrotate(V1horizontalee,-75,'bilinear','crop');
V175ii = imrotate(V1horizontalii,-75,'bilinear','crop');
V1105ee = imrotate(V1horizontalee,-105,'bilinear','crop');
V1105ii = imrotate(V1horizontalii,-105,'bilinear','crop');
V1120ee = imrotate(V1horizontalee,-120,'bilinear','crop');
V1120ii = imrotate(V1horizontalii,-120,'bilinear','crop');
V1135ee = imrotate(V1horizontalee,-135,'bilinear','crop');
V1135ii = imrotate(V1horizontalii,-135,'bilinear','crop');
V1150ee = imrotate(V1horizontalee,-150,'bilinear','crop');
V1150ii = imrotate(V1horizontalii,-150,'bilinear','crop');
V1165ee = imrotate(V1horizontalee,-165,'bilinear','crop');
V1165ii = imrotate(V1horizontalii,-165,'bilinear','crop');


%  n= -1;
%  n15 = -1;
%  n30 = -1;
%  n45 = -1;
%  n60 = -1;
%  n75 = -1;
%  n90 = -1;
%  n105 = -1;
%  n120 = -1;
%  n135 = -1;
%  n150 = -1;
%  n165 = -1;
% 
% for i = 1:V1e  % this is to optimize kernel parameters of receptive field
%   
%     b = sum(V1horizontalii(i,:));
%     
%     
%     b15 = sum(V115ii(i,:));
%     
%     b30 = sum(V130ii(i,:));
%     b45 = sum(V145ii(i,:));
%     b60 = sum(V160ii(i,:));
%     b75 = sum(V175ii(i,:));
%     b90 = sum(V1verticalii(i,:));
%     b105 = sum(V1105ii(i,:));
%     b120 = sum(V1120ii(i,:));
%     b135 = sum(V1135ii(i,:));
%     b150 = sum(V1150ii(i,:));
%     b165 = sum(V1165ii(i,:));
%     if b == 0 
%        
%            n = n + 1;
%         V1horizontali(i-n,:)=[];               
%     end
%     
%     if b15 == 0
%         n15 = n15 + 1;
%         V115i(i-n15,:)=[];
%     end
%     if b30 == 0
%         n30 = n30 + 1;
%         V130i(i-n30,:)=[];
%     end
%     
%     if b45 == 0
%         n45 = n45 + 1;
%         V145i(i-n45,:)=[];
%     end
%     
%     if b60 == 0
%         n60 = n60 + 1;
%         V160i(i-n60,:)=[];
%     end
%     
%     if b75 == 0
%         n75 = n75 + 1;
%         V175i(i-n75,:)=[];
%     end
%     
%     if b90 == 0
%         n90 = n90 + 1;
%         V1verticali(i-n90,:)=[];
%     end
%     
%     if b105 == 0
%         n105 = n105 + 1;
%         V1105i(i-n105,:)=[];
%     end
%     
%     if b120 == 0
%         n120 = n120 + 1;
%         V1120i(i-n120,:)=[];
%     end
%     
%     if b135 == 0
%         n135 = n135 + 1;
%         V1135i(i-n135,:)=[];
%     end
%     
%     if b150 == 0
%         n150 = n150 + 1;
%         V1150i(i-n150,:)=[];
%     end
%     
%     if b165 == 0
%         n165 = n165 + 1;
%         V1165i(i-n165,:)=[];
%     end
% end
% 
% n= -1;
%  n15 = -1;
%  n30 = -1;
%  n45 = -1;
%  n60 = -1;
%  n75 = -1;
%  n90 = -1;
%  n105 = -1;
%  n120 = -1;
%  n135 = -1;
%  n150 = -1;
%  n165 = -1;
% 
% for i = 1:V2e  % this is to optimize kernel parameters of receptive field
%   
%     b = sum(V2horizontalii(i,:));
%     
%     
%     b15 = sum(V215ii(i,:));
%     
%     b30 = sum(V230ii(i,:));
%     b45 = sum(V245ii(i,:));
%     b60 = sum(V260ii(i,:));
%     b75 = sum(V275ii(i,:));
%     b90 = sum(V2verticalii(i,:));
%     b105 = sum(V2105ii(i,:));
%     b120 = sum(V2120ii(i,:));
%     b135 = sum(V2135ii(i,:));
%     b150 = sum(V2150ii(i,:));
%     b165 = sum(V2165ii(i,:));
%     if b == 0 
%        
%            n = n + 1;
%         V2horizontali(i-n,:)=[];               
%     end
%     
%     if b15 == 0
%         n15 = n15 + 1;
%         V215i(i-n15,:)=[];
%     end
%     if b30 == 0
%         n30 = n30 + 1;
%         V230i(i-n30,:)=[];
%     end
%     
%     if b45 == 0
%         n45 = n45 + 1;
%         V245i(i-n45,:)=[];
%     end
%     
%     if b60 == 0
%         n60 = n60 + 1;
%         V260i(i-n60,:)=[];
%     end
%     
%     if b75 == 0
%         n75 = n75 + 1;
%         V275i(i-n75,:)=[];
%     end
%     
%     if b90 == 0
%         n90 = n90 + 1;
%         V2verticali(i-n90,:)=[];
%     end
%     
%     if b105 == 0
%         n105 = n105 + 1;
%         V2105i(i-n105,:)=[];
%     end
%     
%     if b120 == 0
%         n120 = n120 + 1;
%         V2120i(i-n120,:)=[];
%     end
%     
%     if b135 == 0
%         n135 = n135 + 1;
%         V2135i(i-n135,:)=[];
%     end
%     
%     if b150 == 0
%         n150 = n150 + 1;
%         V2150i(i-n150,:)=[];
%     end
%     
%     if b165 == 0
%         n165 = n165 + 1;
%         V2165i(i-n165,:)=[];
%     end
% end
% 
% k= -1;
% k15 = -1;
%  k30 = -1;
%  k45 = -1;
%  k60 = -1;
%  k75 = -1;
%  k90 = -1;
%  k105 = -1;
%  k120 = -1;
%  k135 = -1;
%  k150 = -1;
%  k165 = -1;
% for i = 1:V1e   % this is to optimize kernel parameters of receptive field
%   
%     h = sum(V1horizontalii(:,i));
%     disp(h) 
%     disp(k)
%     
%      h15 = sum(V115ii(:,i));
%     
%     h30 = sum(V130ii(:,i));
%     h45 = sum(V145ii(:,i));
%     h60 = sum(V160ii(:,i));
%     h75 = sum(V175ii(:,i));
%     h90 = sum(V1verticalii(:,i));
%     h105 = sum(V1105ii(:,i));
%     h120 = sum(V1120ii(:,i));
%     h135 = sum(V1135ii(:,i));
%     h150 = sum(V1150ii(:,i));
%     h165 = sum(V1165ii(:,i));
%     
%     if h == 0 
%        
%            k = k + 1;
%         V1horizontali(:,i-k)=[];
%            
%         
%     end
%     
%     if h15 == 0
%         k15 = k15 + 1;
%         V115i(:,i-k15)=[];
%     end
%     if h30 == 0
%         k30 = k30 + 1;
%         V130i(:,i-k30)=[];
%     end
%     
%     if h45 == 0
%         k45 = k45 + 1;
%         V145i(:,i-k45)=[];
%     end
%     
%     if h60 == 0
%         k60 = k60 + 1;
%         V160i(:,i-k60)=[];
%     end
%     
%     if h75 == 0
%         k75 = k75 + 1;
%         V175i(:,i-k75)=[];
%     end
%     
%     if h90 == 0
%         k90 = k90 + 1;
%         V1verticali(:,i-k90)=[];
%     end
%     
%     if h105 == 0
%         k105 = k105 + 1;
%         V1105i(:,i-k105)=[];
%     end
%     
%     if h120 == 0
%         k120 = k120 + 1;
%         V1120i(:,i-k120)=[];
%     end
%     
%     if h135 == 0
%         k135 = k135 + 1;
%         V1135i(:,i-k135)=[];
%     end
%     
%     if h150 == 0
%         k150 = k150 + 1;
%         V1150i(:,i-k150)=[];
%     end
%     
%     if h165 == 0
%         k165 = k165 + 1;
%         V1165i(:,i-k165)=[];
%     end
%     
% end
% 
% k= -1;
% k15 = -1;
%  k30 = -1;
%  k45 = -1;
%  k60 = -1;
%  k75 = -1;
%  k90 = -1;
%  k105 = -1;
%  k120 = -1;
%  k135 = -1;
%  k150 = -1;
%  k165 = -1;
% for i = 1:V2e   % this is to optimize kernel parameters of receptive field
%   
%     h = sum(V2horizontalii(:,i));
%     disp(h) 
%     disp(k)
%     
%      h15 = sum(V215ii(:,i));
%     
%     h30 = sum(V230ii(:,i));
%     h45 = sum(V245ii(:,i));
%     h60 = sum(V260ii(:,i));
%     h75 = sum(V275ii(:,i));
%     h90 = sum(V2verticalii(:,i));
%     h105 = sum(V2105ii(:,i));
%     h120 = sum(V2120ii(:,i));
%     h135 = sum(V2135ii(:,i));
%     h150 = sum(V2150ii(:,i));
%     h165 = sum(V2165ii(:,i));
%     
%     if h == 0 
%            k = k + 1;
%         V2horizontali(:,i-k)=[];  
%     end
%     if h15 == 0
%         k15 = k15 + 1;
%         V215i(:,i-k15)=[];
%     end
%     if h30 == 0
%         k30 = k30 + 1;
%         V230i(:,i-k30)=[];
%     end
%     if h45 == 0
%         k45 = k45 + 1;
%         V245i(:,i-k45)=[];
%     end 
%     if h60 == 0
%         k60 = k60 + 1;
%         V260i(:,i-k60)=[];
%     end
%     if h75 == 0
%         k75 = k75 + 1;
%         V275i(:,i-k75)=[];
%     end
%     if h90 == 0
%         k90 = k90 + 1;
%         V2verticali(:,i-k90)=[];
%     end
%     if h105 == 0
%         k105 = k105 + 1;
%         V2105i(:,i-k105)=[];
%     end
%     if h120 == 0
%         k120 = k120 + 1;
%         V2120i(:,i-k120)=[];
%     end
%     if h135 == 0
%         k135 = k135 + 1;
%         V2135i(:,i-k135)=[];
%     end
%     if h150 == 0
%         k150 = k150 + 1;
%         V2150i(:,i-k150)=[];
%     end
%     if h165 == 0
%         k165 = k165 + 1;
%         V2165i(:,i-k165)=[];
%     end
% end
% 
% z= -1;
%  z15 = -1;
%  z30 = -1;
%  z45 = -1;
%  z60 = -1;
%  z75 = -1;
%  z90 = -1;
%  z105 = -1;
%  z120 = -1;
%  z135 = -1;
%  z150 = -1;
%  z165 = -1;
% for i = 1:V1e  % this is to optimize kernel parameters of receptive field
%   
%     b = sum(V1horizontalee(i,:));
%     disp(b) 
%     disp(i)
%     
%     b15 = sum(V115ee(i,:));
%     
%     b30 = sum(V130ee(i,:));
%     b45 = sum(V145ee(i,:));
%     b60 = sum(V160ee(i,:));
%     b75 = sum(V175ee(i,:));
%     b90 = sum(V1verticalee(i,:));
%     b105 = sum(V1105ee(i,:));
%     b120 = sum(V1120ee(i,:));
%     b135 = sum(V1135ee(i,:));
%     b150 = sum(V1150ee(i,:));
%     b165 = sum(V1165ee(i,:));
%     if b == 0 
%        
%            z = z + 1;
%         V1horizontale(i-z,:)=[];               
%     end
%     
%     if b15 == 0
%         z15 = z15 + 1;
%         V115e(i-z15,:)=[];
%     end
%     if b30 == 0
%         z30 = z30 + 1;
%         V130e(i-z30,:)=[];
%     end
%     
%     if b45 == 0
%         z45 = z45 + 1;
%         V145e(i-z45,:)=[];
%     end
%     
%     if b60 == 0
%         z60 = z60 + 1;
%         V160e(i-z60,:)=[];
%     end
%     
%     if b75 == 0
%         z75 = z75 + 1;
%         V175e(i-z75,:)=[];
%     end
%     
%     if b90 == 0
%         z90 = z90 + 1;
%         V1verticale(i-z90,:)=[];
%     end
%     
%     if b105 == 0
%         z105 = z105 + 1;
%         V1105e(i-z105,:)=[];
%     end
%     
%     if b120 == 0
%         z120 = z120 + 1;
%         V1120e(i-z120,:)=[];
%     end
%     
%     if b135 == 0
%         z135 = z135 + 1;
%         V1135e(i-z135,:)=[];
%     end
%     
%     if b150 == 0
%         z150 = z150 + 1;
%         V1150e(i-z150,:)=[];
%     end
%     
%     if b165 == 0
%         z165 = z165 + 1;
%         V1165e(i-z165,:)=[];
%     end
% end
% 
% z= -1;
%  z15 = -1;
%  z30 = -1;
%  z45 = -1;
%  z60 = -1;
%  z75 = -1;
%  z90 = -1;
%  z105 = -1;
%  z120 = -1;
%  z135 = -1;
%  z150 = -1;
%  z165 = -1;
% for i = 1:V2e  % this is to optimize kernel parameters of receptive field
%   
%     b = sum(V2horizontalee(i,:));
%     disp(b) 
%     disp(i)
%     
%     b15 = sum(V215ee(i,:));
%     b30 = sum(V230ee(i,:));
%     b45 = sum(V245ee(i,:));
%     b60 = sum(V260ee(i,:));
%     b75 = sum(V275ee(i,:));
%     b90 = sum(V2verticalee(i,:));
%     b105 = sum(V2105ee(i,:));
%     b120 = sum(V2120ee(i,:));
%     b135 = sum(V2135ee(i,:));
%     b150 = sum(V2150ee(i,:));
%     b165 = sum(V2165ee(i,:));
%     if b == 0  
%            z = z + 1;
%         V2horizontale(i-z,:)=[];               
%     end
%     if b15 == 0
%         z15 = z15 + 1;
%         V215e(i-z15,:)=[];
%     end
%     if b30 == 0
%         z30 = z30 + 1;
%         V230e(i-z30,:)=[];
%     end
%     if b45 == 0
%         z45 = z45 + 1;
%         V245e(i-z45,:)=[];
%     end
%     if b60 == 0
%         z60 = z60 + 1;
%         V260e(i-z60,:)=[];
%     end
%     if b75 == 0
%         z75 = z75 + 1;
%         V275e(i-z75,:)=[];
%     end
%     if b90 == 0
%         z90 = z90 + 1;
%         V2verticale(i-z90,:)=[];
%     end
%     if b105 == 0
%         z105 = z105 + 1;
%         V2105e(i-z105,:)=[];
%     end
%     if b120 == 0
%         z120 = z120 + 1;
%         V2120e(i-z120,:)=[];
%     end
%     if b135 == 0
%         z135 = z135 + 1;
%         V2135e(i-z135,:)=[];
%     end
%     if b150 == 0
%         z150 = z150 + 1;
%         V2150e(i-z150,:)=[];
%     end
%     if b165 == 0
%         z165 = z165 + 1;
%         V2165e(i-z165,:)=[];
%     end
% end
% 
% o= -1;
% o15 = -1;
%  o30 = -1;
%  o45 = -1;
%  o60 = -1;
%  o75 = -1;
%  o90 = -1;
%  o105 = -1;
%  o120 = -1;
%  o135 = -1;
%  o150 = -1;
%  o165 = -1;
% for i = 1:V1e   % this is to optimize kernel parameters of receptive field
%   
%     h = sum(V1horizontalee(:,i));
%     disp(h) 
%     
%     
%      h15 = sum(V115ee(:,i));
%     
%     h30 = sum(V130ee(:,i));
%     h45 = sum(V145ee(:,i));
%     h60 = sum(V160ee(:,i));
%     h75 = sum(V175ee(:,i));
%     h90 = sum(V1verticalee(:,i));
%     h105 = sum(V1105ee(:,i));
%     h120 = sum(V1120ee(:,i));
%     h135 = sum(V1135ee(:,i));
%     h150 = sum(V1150ee(:,i));
%     h165 = sum(V1165ee(:,i));
%     
%     if h == 0 
%        
%            o = o + 1;
%         V1horizontale(:,i-o)=[];
%            
%         
%     end
%     
%     if h15 == 0
%         o15 = o15 + 1;
%         V115e(:,i-o15)=[];
%     end
%     if h30 == 0
%         o30 = o30 + 1;
%         V130e(:,i-o30)=[];
%     end
%     
%     if h45 == 0
%         o45 = o45 + 1;
%         V145e(:,i-o45)=[];
%     end
%     
%     if h60 == 0
%         o60 = o60 + 1;
%         V160e(:,i-o60)=[];
%     end
%     
%     if h75 == 0
%         o75 = o75 + 1;
%         V175e(:,i-o75)=[];
%     end
%     
%     if h90 == 0
%         o90 = o90 + 1;
%         V1verticale(:,i-o90)=[];
%     end
%     
%     if h105 == 0
%         o105 = o105 + 1;
%         V1105e(:,i-o105)=[];
%     end
%     
%     if h120 == 0
%         o120 = o120 + 1;
%         V1120e(:,i-o120)=[];
%     end
%     
%     if h135 == 0
%         o135 = o135 + 1;
%         V1135e(:,i-o135)=[];
%     end
%     
%     if h150 == 0
%         o150 = o150 + 1;
%         V1150e(:,i-o150)=[];
%     end
%     
%     if h165 == 0
%         o165 = o165 + 1;
%         V1165e(:,i-o165)=[];
%     end
%     
% end
% 
% o= -1;
% o15 = -1;
%  o30 = -1;
%  o45 = -1;
%  o60 = -1;
%  o75 = -1;
%  o90 = -1;
%  o105 = -1;
%  o120 = -1;
%  o135 = -1;
%  o150 = -1;
%  o165 = -1;
% for i = 1:V2e   % this is to optimize kernel parameters of receptive field
%     h = sum(V2horizontalee(:,i));
%     disp(h) 
%      h15 = sum(V215ee(:,i));
%     h30 = sum(V230ee(:,i));
%     h45 = sum(V245ee(:,i));
%     h60 = sum(V260ee(:,i));
%     h75 = sum(V275ee(:,i));
%     h90 = sum(V2verticalee(:,i));
%     h105 = sum(V2105ee(:,i));
%     h120 = sum(V2120ee(:,i));
%     h135 = sum(V2135ee(:,i));
%     h150 = sum(V2150ee(:,i));
%     h165 = sum(V2165ee(:,i));
%     if h == 0 
%            o = o + 1;
%         V2horizontale(:,i-o)=[];
%     end
%     if h15 == 0
%         o15 = o15 + 1;
%         V215e(:,i-o15)=[];
%     end
%     if h30 == 0
%         o30 = o30 + 1;
%         V230e(:,i-o30)=[];
%     end
%     if h45 == 0
%         o45 = o45 + 1;
%         V245e(:,i-o45)=[];
%     end
%     if h60 == 0
%         o60 = o60 + 1;
%         V260e(:,i-o60)=[];
%     end
%     if h75 == 0
%         o75 = o75 + 1;
%         V275e(:,i-o75)=[];
%     end
%     if h90 == 0
%         o90 = o90 + 1;
%         V2verticale(:,i-o90)=[];
%     end
%     if h105 == 0
%         o105 = o105 + 1;
%         V2105e(:,i-o105)=[];
%     end
%     if h120 == 0
%         o120 = o120 + 1;
%         V2120e(:,i-o120)=[];
%     end
%     if h135 == 0
%         o135 = o135 + 1;
%         V2135e(:,i-o135)=[];
%     end
%     if h150 == 0
%         o150 = o150 + 1;
%         V2150e(:,i-o150)=[];
%     end
%     if h165 == 0
%         o165 = o165 + 1;
%         V2165e(:,i-o165)=[];
%     end
% end



% ellipticalfield = V1horizontale-V1horizontali;

% ellipticalfieldv = V1verticale - V1verticali;
% concentricV1e=fspecial('gaussian',[V1e, V1e], 4); % excitation normal gaussian 
% concentricV1i=fspecial('gaussian',[V1e, V1e], 8); % inhibition 

Weights_xy3_e=fspecial('gaussian',[V1e V1e], 4); % excitation normal gaussian 
Weights_xy3_i=fspecial('gaussian',[V1e V1e], 8);
%of inputs to neurons in question
%need to evaluate new inputs based on synaptic weights later
% OFFRGCe=fspecial('gaussian',[25 25], 1 ); % excitation normal gaussian 
% OFFRGCi=fspecial('gaussian',[25 25], 2); % inhibition 
% 
% OFFLGNe=fspecial('gaussian',[35 35], 6); % excitation normal gaussian 
% OFFLGNi=fspecial('gaussian',[35 35], 3); % inhibition 
% 
% OFFV1e=fspecial('gaussian',[54 54], 8); % excitation normal gaussian 
% OFFV1i=fspecial('gaussian',[54 54], 4); % inhibition 
% 
% OFFfeedbacke=fspecial('gaussian',[54 54], 8); % excitation normal gaussian 
% OFFfeedbacki=fspecial('gaussian',[54 54], 4);

% xa=0.55; 
xa = 0.55;
   %leak level/decay contsants for x equation

% Iex = shiftdim(Ie,1);
% Iixx = shiftdim(Iix,1);
 

 RGC=(Ie);
 RGC1=RGC;
 RGC(:,:,:) = 0;
 RGC1(:,:,:)=0;
 LGN=RGC;
 LGN(:,:,:)=0;
 LGN1=LGN;
 LGN1(:,:,:)=0;
 LGN2=LGN1;
 LGN2(:,:,:)=0;
 RGC2=RGC1;
 RGC2(:,:,:)=0;


V1 = LGN*0;
V1x = zeros(size(V1,1),size(V1,2));
V1y =zeros(size(V1,1),size(V1,2));
V115x=zeros(size(V1,1),size(V1,2));
V115y= zeros(size(V1,1),size(V1,2));
V130x= zeros(size(V1,1),size(V1,2));
V130y= zeros(size(V1,1),size(V1,2));
V145x= zeros(size(V1,1),size(V1,2));
V145y= zeros(size(V1,1),size(V1,2));
V160x= zeros(size(V1,1),size(V1,2));
V160y= zeros(size(V1,1),size(V1,2));
V175x= zeros(size(V1,1),size(V1,2));
V175y= zeros(size(V1,1),size(V1,2));
V190x= zeros(size(V1,1),size(V1,2));
V190y= zeros(size(V1,1),size(V1,2));
V1105x= zeros(size(V1,1),size(V1,2));
V1105y= zeros(size(V1,1),size(V1,2));
V1120x= zeros(size(V1,1),size(V1,2));
V1120y= zeros(size(V1,1),size(V1,2));
V1135x= zeros(size(V1,1),size(V1,2));
V1135y= zeros(size(V1,1),size(V1,2));
V1150x= zeros(size(V1,1),size(V1,2));
V1150y= zeros(size(V1,1),size(V1,2));
V1165x= zeros(size(V1,1),size(V1,2));
V1165y= zeros(size(V1,1),size(V1,2));
V15 = V1*0;
V110 = V1*0;
V120 = V1*0;
V125 = V1*0;
V135 = V1*0;
V140 = V1*0;
V150 = V1*0;
V155 = V1*0;
V165 = V1*0;
V170 = V1*0;
V180 = V1*0;
V185 = V1*0;
V195 = V1*0;
V1100 = V1*0;
V1110 = V1*0;
V1115 = V1*0;
V1125 = V1*0;
V1130 = V1*0;
V1140 = V1*0;
V1145 = V1*0;
V1155 = V1*0;
V1160 = V1*0;
V1170 = V1*0;
V1175 = V1*0;

V25 = V1*0;
V210 = V1*0;
V220 = V1*0;
V225 = V1*0;
V235 = V1*0;
V240 = V1*0;
V250 = V1*0;
V255 = V1*0;
V265 = V1*0;
V270 = V1*0;
V280 = V1*0;
V285 = V1*0;
V295 = V1*0;
V2100 = V1*0;
V2110 = V1*0;
V2115 = V1*0;
V2125 = V1*0;
V2130 = V1*0;
V2140 = V1*0;
V2145 = V1*0;
V2155 = V1*0;
V2160 = V1*0;
V2170 = V1*0;
V2175 = V1*0;

V145 = V1;
V145(:,:,:)=0;
V1135 = V1;
V1135(:,:,:) = 0;
V1(:,:,:)=0;
V190 = V1;
V190(:,:,:)=0;
V115 = V1;
V115(:,:,:) = 0;
V130 = V1;
V130(:,:,:) = 0;
V160=V1;
V160(:,:,:)=0;
V175=V1;
V175(:,:,:)=0;
V1105=V1;
V1105(:,:,:)=0;
V1120=V1;
V1120(:,:,:)=0;
V1150=V1;
V1150(:,:,:)=0;
V1165=V1;
V1165(:,:,:)=0;

V2 = LGN;
V245 = V1;
V245(:,:,:)=0;
V2135 = V1;
V2135(:,:,:) = 0;
V2(:,:,:)=0;
V290 = V1;
V290(:,:,:)=0;
V215 = V1;
V215(:,:,:) = 0;
V230 = V1;
V230(:,:,:) = 0;
V260=V1;
V260(:,:,:)=0;
V275=V1;
V275(:,:,:)=0;
V2105=V1;
V2105(:,:,:)=0;
V2120=V1;
V2120(:,:,:)=0;
V2150=V1;
V2150(:,:,:)=0;
V2165=V1;
V2165(:,:,:)=0;

V2090=V1;
V2090(:,:,:)=0;
V215105=V1;
V215105(:,:,:)=0;
V230120=V1;
V230120(:,:,:)=0;
V245135=V1;
V245135(:,:,:)=0;
V260150=V1;
V260150(:,:,:)=0;
V275165=V1;
V275165(:,:,:)=0;
V1feedback=V1;
V1feedback(:,:,:)=0;


RGCnumwidth = size(RGC,1);
RGCnumlength = RGCnumwidth;

fs = 1;% feedback strength
is = 1; % interaction strength

% for i= 1:n     
%    
%     dx= -xa.*RGC(:,:,i)+...
%         1.*(1-RGC(:,:,i)).*(convn(Ie(:,:,i), Weights_xy_e, 'same'))-...
%         1.*(1+RGC(:,:,i)).*(convn(Iix(:,:,i), Weights_xy_i, 'same'));
%       
%    RGC(:,:,i+1)=RGC(:,:,i)+dT.*dx;
%    
%    dx1= -xa.*LGN(:,:,i)+...
%         1.*(1-LGN(:,:,i)).*(convn(RGC(:,:,i), Weights_xy1_e, 'same')+fs.*(convn(V1(:,:,i),Weights_xy3_e,'same')))-...
%         1.*(1+LGN(:,:,i)).*(convn(RGC(:,:,i)./1.16, Weights_xy1_i, 'same')-fs.*(convn(V1(:,:,i)./1.16,Weights_xy3_i,'same')));
%       
%    LGN(:,:,i+1)=LGN(:,:,i)+dT.*dx1;
%    
%    dx2= -xa.*V1(:,:,i)+...
%         1.*(1-V1(:,:,i)).*(convn(LGN(:,:,i), Weights_xy2_e, 'same'))-...
%         1.*(1+V1(:,:,i)).*(convn(LGN(:,:,i)./1.16, Weights_xy2_i, 'same'));
%       
%    V1(:,:,i+1)=V1(:,:,i)+dT.*dx2;
%      
%      % equation to model neuronal activity based on previous responses
% end
 squwidth = 1;
    disp(squwidth)
    squlength = squwidth;
    numwidth = floor(size(V1,2)/squlength);
    numlength = floor(size(V1,1)/squwidth);
    columnbox = zeros(numwidth,12);
    rowbox = zeros(12,numlength);
    meanbox1 = zeros(numlength,numwidth,36);
    meanbox2 = zeros(numlength,numwidth,12);
    comparebox = zeros(numlength,numwidth,12);
    negabox = zeros(numlength,numwidth,12);
 Verticalscanbox = meanbox1;
 Verticalscanbox(:,:,:) = 0;
 Horizontalscanbox = Verticalscanbox*0;
 localinteractionbox = zeros(numlength,numwidth,72,time);
 localinteraction = zeros(numlength,numwidth,72,time);
 
 e=(fspecial('gaussian',[1,12], 1)); % excitation normal gaussian 
ii=(fspecial('gaussian',[1,12], 1.6));

 ee=(fspecial('gaussian',[12,12], 1)); % excitation normal gaussian 
iii=(fspecial('gaussian',[12,12], 1.6));

inhibconstant = 1;
inhibconstanti = 1;

Interactioni = Interactioni*inhibconstanti;
RGCii = RGCii*inhibconstant;
LGNii = LGNii*inhibconstant;
V1verticali=V1verticali*inhibconstant;
V1horizontali=V1horizontali*inhibconstant;

V110i = V110i*inhibconstant;
V120i = V120i*inhibconstant;
V130i = V130i*inhibconstant;
V140i = V140i*inhibconstant;
V150i = V150i*inhibconstant;
V160i = V160i*inhibconstant;
V170i = V170i*inhibconstant;
V180i = V180i*inhibconstant;
V1100i = V1100i*inhibconstant;
V1110i = V1110i*inhibconstant;
V1120i = V1120i*inhibconstant;
V1130i = V1130i*inhibconstant;
V1140i = V1140i*inhibconstant;
V1150i = V1150i*inhibconstant;
V1160i = V1160i*inhibconstant;
V1170i = V1170i*inhibconstant;
V2verticali=V2verticali*inhibconstant;
V2horizontali=V2horizontali*inhibconstant;
V210i = V210i*inhibconstant;
V220i = V220i*inhibconstant;
V230i = V230i*inhibconstant;
V240i = V240i*inhibconstant;
V250i = V250i*inhibconstant;
V260i = V260i*inhibconstant;
V270i = V270i*inhibconstant;
V280i = V280i*inhibconstant;
V2100i = V2100i*inhibconstant;
V2110i = V2110i*inhibconstant;
V2120i = V2120i*inhibconstant;
V2130i = V2130i*inhibconstant;
V2140i = V2140i*inhibconstant;
V2150i = V2150i*inhibconstant;
V2160i = V2160i*inhibconstant;
V2170i = V2170i*inhibconstant;

exciteconstant = 1;
exciteconstanti = 1;
Interactione = Interactione*exciteconstanti;
RGCee = RGCee*exciteconstant;
LGNee = LGNee*exciteconstant;
V1verticale=V1verticale*exciteconstant;
V1horizontale=V1horizontale*exciteconstant;
V110e = V110e*exciteconstant;
V120e = V120e*exciteconstant;
V130e = V130e*exciteconstant;
V140e = V140e*exciteconstant;
V150e = V150e*exciteconstant;
V160e = V160e*exciteconstant;
V170e = V170e*exciteconstant;
V180e = V180e*exciteconstant;
V1100e = V1100e*exciteconstant;
V1110e = V1110e*exciteconstant;
V1120e = V1120e*exciteconstant;
V1130e = V1130e*exciteconstant;
V1140e = V1140e*exciteconstant;
V1150e = V1150e*exciteconstant;
V1160e = V1160e*exciteconstant;
V1170e = V1170e*exciteconstant;
V2verticale=V2verticale*exciteconstant;
V2horizontale=V2horizontale*exciteconstant;
V210e = V210e*exciteconstant;
V220e = V220e*exciteconstant;
V230e = V230e*exciteconstant;
V240e = V240e*exciteconstant;
V250e = V250e*exciteconstant;
V260e = V260e*exciteconstant;
V270e = V270e*exciteconstant;
V280e = V280e*exciteconstant;
V2100e = V2100e*exciteconstant;
V2110e = V2110e*exciteconstant;
V2120e = V2120e*exciteconstant;
V2130e = V2130e*exciteconstant;
V2140e = V2140e*exciteconstant;
V2150e = V2150e*exciteconstant;
V2160e = V2160e*exciteconstant;
V2170e = V2170e*exciteconstant;

for i= 1:time     %
    
        disp(i)
        
%         dRGC= -xa.*RGC(:,:,i)+...
%         1.*(1-RGC(:,:,i)).*(convn(Ie(:,:,i), RGCee, 'same'))-...
%         1.*(1+RGC(:,:,i)).*(convn(Ie(:,:,i), RGCii, 'same'));
%       
%    RGC(:,:,i+1)=RGC(:,:,i)+dT.*dRGC;
% 
% dLGN1= -xa.*LGN(:,:,i)+...
%          1.*(1-LGN(:,:,i)).*(convn(RGC(:,:,i), LGNee, 'same'))-...
%          1.*(1+LGN(:,:,i)).*(convn(RGC(:,:,i), LGNii, 'same'));
%        
%    LGN(:,:,i+1)=LGN(:,:,i)+dT.*dLGN1;

   dV190= -xa.*V190(:,:,i)+...
        1.*(1-V190(:,:,i)).*fs.*(convn(Ie(:,:,i)+V290(:,:,i), V1verticale, 'same'))-...
        1.*(1+V190(:,:,i)).*fs.*(convn(Ie(:,:,i)+V290(:,:,i), V1verticali, 'same'));
   % V190(:,:,i+1)=V190(:,:,i)+dT.*dV190;
    
    dV1= -xa.*V1(:,:,i)+...
        1.*(1-V1(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2(:,:,i), V1horizontale, 'same'))-...
        1.*(1+V1(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2(:,:,i), V1horizontali, 'same'));
   % V1(:,:,i+1)=V1(:,:,i)+dT.*dV1;
   
      
           dV110= -xa.*V110(:,:,i)+...
        1.*(1-V110(:,:,i)).*fs.*(convn(Ie(:,:,i)+V210(:,:,i), V110e, 'same'))-...
        1.*(1+V110(:,:,i)).*fs.*(convn(Ie(:,:,i)+V210(:,:,i), V110i, 'same'));
    
               dV120= -xa.*V120(:,:,i)+...
        1.*(1-V120(:,:,i)).*fs.*(convn(Ie(:,:,i)+V220(:,:,i), V120e, 'same'))-...
        1.*(1+V120(:,:,i)).*fs.*(convn(Ie(:,:,i)+V220(:,:,i), V120i, 'same'));
    
         
    
                  
                           dV140= -xa.*V140(:,:,i)+...
        1.*(1-V140(:,:,i)).*fs.*(convn(Ie(:,:,i)+V240(:,:,i), V140e, 'same'))-...
        1.*(1+V140(:,:,i)).*fs.*(convn(Ie(:,:,i)+V240(:,:,i), V140i, 'same'));
 
                               dV150= -xa.*V150(:,:,i)+...
        1.*(1-V150(:,:,i)).*fs.*(convn(Ie(:,:,i)+V250(:,:,i), V150e, 'same'))-...
        1.*(1+V150(:,:,i)).*fs.*(convn(Ie(:,:,i)+V250(:,:,i), V150i, 'same'));
    
                                  
    
                                        dV170= -xa.*V170(:,:,i)+...
        1.*(1-V170(:,:,i)).*fs.*(convn(Ie(:,:,i)+V270(:,:,i), V170e, 'same'))-...
        1.*(1+V170(:,:,i)).*fs.*(convn(Ie(:,:,i)+V270(:,:,i), V170i, 'same'));
 
                                         dV180= -xa.*V180(:,:,i)+...
        1.*(1-V180(:,:,i)).*fs.*(convn(Ie(:,:,i)+V280(:,:,i), V180e, 'same'))-...
        1.*(1+V180(:,:,i)).*fs.*(convn(Ie(:,:,i)+V280(:,:,i), V180i, 'same'));
    
                                            
  
                               
    
                                     dV1100= -xa.*V1100(:,:,i)+...
        1.*(1-V1100(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2100(:,:,i), V1100e, 'same'))-...
        1.*(1+V1100(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2100(:,:,i), V1100i, 'same'));
    
                         dV1110= -xa.*V1110(:,:,i)+...
        1.*(1-V1110(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2110(:,:,i), V1110e, 'same'))-...
        1.*(1+V1110(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2110(:,:,i), V1110i, 'same'));
    
                            
    
                                     dV1130= -xa.*V1130(:,:,i)+...
        1.*(1-V1130(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2130(:,:,i), V1130e, 'same'))-...
        1.*(1+V1130(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2130(:,:,i), V1130i, 'same'));
  
                                         dV1140= -xa.*V1140(:,:,i)+...
        1.*(1-V1140(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2140(:,:,i), V1140e, 'same'))-...
        1.*(1+V1140(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2140(:,:,i), V1140i, 'same'));
    
    
    
                                               
                                       dV1160= -xa.*V1160(:,:,i)+...
        1.*(1-V1160(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2160(:,:,i), V1160e, 'same'))-...
        1.*(1+V1160(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2160(:,:,i), V1160i, 'same'));
    
                                       dV1170= -xa.*V1170(:,:,i)+...
        1.*(1-V1170(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2170(:,:,i), V1170e, 'same'))-...
        1.*(1+V1170(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2170(:,:,i), V1170i, 'same'));
    
    
    
   
      
   
    dV130= -xa.*V130(:,:,i)+...
        1.*(1-V130(:,:,i)).*fs.*(convn(Ie(:,:,i)+V230(:,:,i), V130e, 'same'))-...
        1.*(1+V130(:,:,i)).*fs.*(convn(Ie(:,:,i)+V230(:,:,i), V130i, 'same'));
    % V1112d5(:,:,i+1)=V1112d5(:,:,i)+dT.*dV1112d5;
    
    dV160= -xa.*V160(:,:,i)+...
        1.*(1-V160(:,:,i)).*fs.*(convn(Ie(:,:,i)+V260(:,:,i), V160e, 'same'))-...
        1.*(1+V160(:,:,i)).*fs.*(convn(Ie(:,:,i)+V260(:,:,i), V160i, 'same'));
    % V167d5(:,:,i+1)=V167d5(:,:,i)+dT.*dV167d5;
    
    
    
    
    dV1120= -xa.*V1120(:,:,i)+...
        1.*(1-V1120(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2120(:,:,i), V1120e, 'same'))-...
        1.*(1+V1120(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2120(:,:,i), V1120i, 'same'));
    
    dV1150= -xa.*V1150(:,:,i)+...
        1.*(1-V1150(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2150(:,:,i), V1150e, 'same'))-...
        1.*(1+V1150(:,:,i)).*fs.*(convn(Ie(:,:,i)+V2150(:,:,i), V1150i, 'same'));
    
    
    

%    dV190= -xa.*V190(:,:,i)+...
%         1.*(1-V190(:,:,i)).*(convn(LGN(:,:,i)+V290(:,:,i), V1verticale, 'same'))-...
%         1.*(1+V190(:,:,i)).*(convn(LGN(:,:,i)+V290(:,:,i), V1verticali, 'same'));
%    % V190(:,:,i+1)=V190(:,:,i)+dT.*dV190;
%     
%     dV1= -xa.*V1(:,:,i)+...
%         1.*(1-V1(:,:,i)).*(convn(LGN(:,:,i)+V2(:,:,i), V1horizontale, 'same'))-...
%         1.*(1+V1(:,:,i)).*(convn(LGN(:,:,i)+V2(:,:,i), V1horizontali, 'same'));
%    % V1(:,:,i+1)=V1(:,:,i)+dT.*dV1;
%     
%     dV145= -xa.*V145(:,:,i)+...
%         1.*(1-V145(:,:,i)).*(convn(LGN(:,:,i)+V245(:,:,i), V145e, 'same'))-...
%         1.*(1+V145(:,:,i)).*(convn(LGN(:,:,i)+V245(:,:,i), V145i, 'same'));
%    % V145(:,:,i+1)=V145(:,:,i)+dT.*dV145;
%     
%     dV1135= -xa.*V1135(:,:,i)+...
%         1.*(1-V1135(:,:,i)).*(convn(LGN(:,:,i)+V2135(:,:,i), V1135e, 'same'))-...
%         1.*(1+V1135(:,:,i)).*(convn(LGN(:,:,i)+V2135(:,:,i), V1135i, 'same'));
%     % V1135(:,:,i+1)=V1135(:,:,i)+dT.*dV1135;
%       
%     dV115= -xa.*V115(:,:,i)+...
%         1.*(1-V115(:,:,i)).*(convn(LGN(:,:,i)+V215(:,:,i), V115e, 'same'))-...
%         1.*(1+V115(:,:,i)).*(convn(LGN(:,:,i)+V215(:,:,i), V115i, 'same'));
%     % V122d5(:,:,i+1)=V122d5(:,:,i)+dT.*dV122d5;
%     
%     dV130= -xa.*V130(:,:,i)+...
%         1.*(1-V130(:,:,i)).*(convn(LGN(:,:,i)+V230(:,:,i), V130e, 'same'))-...
%         1.*(1+V130(:,:,i)).*(convn(LGN(:,:,i)+V230(:,:,i), V130i, 'same'));
%     % V1112d5(:,:,i+1)=V1112d5(:,:,i)+dT.*dV1112d5;
%     
%     dV160= -xa.*V160(:,:,i)+...
%         1.*(1-V160(:,:,i)).*(convn(LGN(:,:,i)+V260(:,:,i), V160e, 'same'))-...
%         1.*(1+V160(:,:,i)).*(convn(LGN(:,:,i)+V260(:,:,i), V160i, 'same'));
%     % V167d5(:,:,i+1)=V167d5(:,:,i)+dT.*dV167d5;
%     
%     dV175= -xa.*V175(:,:,i)+...
%         1.*(1-V175(:,:,i)).*(convn(LGN(:,:,i)+V275(:,:,i), V175e, 'same'))-...
%         1.*(1+V175(:,:,i)).*(convn(LGN(:,:,i)+V275(:,:,i), V175i, 'same'));
%     % V1157d5(:,:,i+1)=V1157d5(:,:,i)+dT.*dV1157d5;
%     
%     dV1105= -xa.*V1105(:,:,i)+...
%         1.*(1-V1105(:,:,i)).*(convn(LGN(:,:,i)+V2105(:,:,i), V1105e, 'same'))-...
%         1.*(1+V1105(:,:,i)).*(convn(LGN(:,:,i)+V2105(:,:,i), V1105i, 'same'));
%     
%     dV1120= -xa.*V1120(:,:,i)+...
%         1.*(1-V1120(:,:,i)).*(convn(LGN(:,:,i)+V2120(:,:,i), V1120e, 'same'))-...
%         1.*(1+V1120(:,:,i)).*(convn(LGN(:,:,i)+V2120(:,:,i), V1120i, 'same'));
%     
%     dV1150= -xa.*V1150(:,:,i)+...
%         1.*(1-V1150(:,:,i)).*(convn(LGN(:,:,i)+V2150(:,:,i), V1150e, 'same'))-...
%         1.*(1+V1150(:,:,i)).*(convn(LGN(:,:,i)+V2150(:,:,i), V1150i, 'same'));
%     
%     dV1165= -xa.*V1165(:,:,i)+...
%         1.*(1-V1165(:,:,i)).*(convn(LGN(:,:,i)+V2165(:,:,i), V1165e, 'same'))-...
%         1.*(1+V1165(:,:,i)).*(convn(LGN(:,:,i)+V2165(:,:,i), V1165i, 'same'));

%  dV190= -xa.*V290(:,:,i)+...
%         1.*(1-V290(:,:,i)).*(convn(Ie(:,:,i), V1verticale, 'same'))-...
%         1.*(1+V290(:,:,i)).*(convn(Ie(:,:,i), V1verticali, 'same'));
%    % V190(:,:,i+1)=V190(:,:,i)+dT.*dV190;
%     
%     dV1= -xa.*V2(:,:,i)+...
%         1.*(1-V2(:,:,i)).*(convn(Ie(:,:,i), V1horizontale, 'same'))-...
%         1.*(1+V2(:,:,i)).*(convn(Ie(:,:,i), V1horizontali, 'same'));
%    % V1(:,:,i+1)=V1(:,:,i)+dT.*dV1;
%     
%     dV145= -xa.*V245(:,:,i)+...
%         1.*(1-V245(:,:,i)).*(convn(Ie(:,:,i), V145e, 'same'))-...
%         1.*(1+V245(:,:,i)).*(convn(Ie(:,:,i), V145i, 'same'));
%    % V145(:,:,i+1)=V145(:,:,i)+dT.*dV145;
%     
%     dV1135= -xa.*V2135(:,:,i)+...
%         1.*(1-V2135(:,:,i)).*(convn(Ie(:,:,i), V1135e, 'same'))-...
%         1.*(1+V2135(:,:,i)).*(convn(Ie(:,:,i), V1135i, 'same'));
%     % V1135(:,:,i+1)=V1135(:,:,i)+dT.*dV1135;
%       
%     dV115= -xa.*V215(:,:,i)+...
%         1.*(1-V215(:,:,i)).*(convn(Ie(:,:,i), V115e, 'same'))-...
%         1.*(1+V215(:,:,i)).*(convn(Ie(:,:,i), V115i, 'same'));
%     % V122d5(:,:,i+1)=V122d5(:,:,i)+dT.*dV122d5;
%     
%     dV130= -xa.*V230(:,:,i)+...
%         1.*(1-V230(:,:,i)).*(convn(Ie(:,:,i), V130e, 'same'))-...
%         1.*(1+V230(:,:,i)).*(convn(Ie(:,:,i), V130i, 'same'));
%     % V1112d5(:,:,i+1)=V1112d5(:,:,i)+dT.*dV1112d5;
%     
%     dV160= -xa.*V260(:,:,i)+...
%         1.*(1-V260(:,:,i)).*(convn(Ie(:,:,i), V160e, 'same'))-...
%         1.*(1+V260(:,:,i)).*(convn(Ie(:,:,i), V160i, 'same'));
%     % V167d5(:,:,i+1)=V167d5(:,:,i)+dT.*dV167d5;
%     
%     dV175= -xa.*V275(:,:,i)+...
%         1.*(1-V275(:,:,i)).*(convn(Ie(:,:,i), V175e, 'same'))-...
%         1.*(1+V275(:,:,i)).*(convn(Ie(:,:,i), V175i, 'same'));
%     % V1157d5(:,:,i+1)=V1157d5(:,:,i)+dT.*dV1157d5;
%     
%     dV1105= -xa.*V2105(:,:,i)+...
%         1.*(1-V2105(:,:,i)).*(convn(Ie(:,:,i), V1105e, 'same'))-...
%         1.*(1+V2105(:,:,i)).*(convn(Ie(:,:,i), V1105i, 'same'));
%     
%     dV1120= -xa.*V2120(:,:,i)+...
%         1.*(1-V2120(:,:,i)).*(convn(Ie(:,:,i), V1120e, 'same'))-...
%         1.*(1+V2120(:,:,i)).*(convn(Ie(:,:,i), V1120i, 'same'));
%     
%     dV1150= -xa.*V2150(:,:,i)+...
%         1.*(1-V2150(:,:,i)).*(convn(Ie(:,:,i), V1150e, 'same'))-...
%         1.*(1+V2150(:,:,i)).*(convn(Ie(:,:,i), V1150i, 'same'));
%     
%     dV1165= -xa.*V2165(:,:,i)+...
%         1.*(1-V2165(:,:,i)).*(convn(Ie(:,:,i), V1165e, 'same'))-...
%         1.*(1+V2165(:,:,i)).*(convn(Ie(:,:,i), V1165i, 'same'));
    
        V1(:,:,i+1)=V1(:,:,i)+dT.*dV1;
       
        V110(:,:,i+1)=V110(:,:,i)+dT.*dV110;
    
    V120(:,:,i+1)=V120(:,:,i)+dT.*dV120;
    
    V130(:,:,i+1)=V130(:,:,i)+dT.*dV130;
   
    V140(:,:,i+1)=V140(:,:,i)+dT.*dV140;
    
    V150(:,:,i+1)=V150(:,:,i)+dT.*dV150;
    
    V160(:,:,i+1)=V160(:,:,i)+dT.*dV160;
   
    V170(:,:,i+1)=V170(:,:,i)+dT.*dV170;
    
    V180(:,:,i+1)=V180(:,:,i)+dT.*dV180;
    
    V190(:,:,i+1)=V190(:,:,i)+dT.*dV190;
    
    V1100(:,:,i+1)=V1100(:,:,i)+dT.*dV1100;
    
    V1110(:,:,i+1)=V1110(:,:,i)+dT.*dV1110;
    
    V1120(:,:,i+1)=V1120(:,:,i)+dT.*dV1120;
    
    V1130(:,:,i+1)=V1130(:,:,i)+dT.*dV1130;
   
    V1140(:,:,i+1)=V1140(:,:,i)+dT.*dV1140;
    
    V1150(:,:,i+1)=V1150(:,:,i)+dT.*dV1150;
   
    V1160(:,:,i+1)=V1160(:,:,i)+dT.*dV1160;
   
    V1170(:,:,i+1)=V1170(:,:,i)+dT.*dV1170;
   
    
%     V1(:,:,i+1)=V1(:,:,i)+dT.*dV1;
%     V115(:,:,i+1)=V115(:,:,i)+dT.*dV115;
%     V130(:,:,i+1)=V130(:,:,i)+dT.*dV130;
%     V145(:,:,i+1)=V145(:,:,i)+dT.*dV145;
%     V160(:,:,i+1)=V160(:,:,i)+dT.*dV160;
%     V175(:,:,i+1)=V175(:,:,i)+dT.*dV175;
%     V190(:,:,i+1)=V190(:,:,i)+dT.*dV190;
%     V1105(:,:,i+1)=V1105(:,:,i)+dT.*dV1105;
%     V1120(:,:,i+1)=V1120(:,:,i)+dT.*dV1120;
%     V1135(:,:,i+1)=V1135(:,:,i)+dT.*dV1135;
%     V1150(:,:,i+1)=V1150(:,:,i)+dT.*dV1150;
%     V1165(:,:,i+1)=V1165(:,:,i)+dT.*dV1165;
    
%     localinteractionbox(:,:,1,i) = V1(:,:,i)+dT.*dV1; % % four dimensions because three dimensional matrix + time
%     localinteractionbox(:,:,2,i) = V115(:,:,i)+dT.*dV115;
%     localinteractionbox(:,:,3,i) = V130(:,:,i)+dT.*dV130;
%     localinteractionbox(:,:,4,i) = V145(:,:,i)+dT.*dV145;
%     localinteractionbox(:,:,5,i) = V160(:,:,i)+dT.*dV160;
%     localinteractionbox(:,:,6,i) = V175(:,:,i)+dT.*dV175;
%     localinteractionbox(:,:,7,i) = V190(:,:,i)+dT.*dV190;
%     localinteractionbox(:,:,8,i) = V1105(:,:,i)+dT.*dV1105;
%     localinteractionbox(:,:,9,i) = V1120(:,:,i)+dT.*dV1120;
%     localinteractionbox(:,:,10,i) = V1135(:,:,i)+dT.*dV1135;
%     localinteractionbox(:,:,11,i) = V1150(:,:,i)+dT.*dV1150;
%     localinteractionbox(:,:,12,i) = V1165(:,:,i)+dT.*dV1165;
%     localinteractionbox(:,:,13,i) = V1(:,:,i)+dT.*dV1;
%     localinteractionbox(:,:,14,i) = V115(:,:,i)+dT.*dV115;
%     localinteractionbox(:,:,15,i) = V130(:,:,i)+dT.*dV130;
%     localinteractionbox(:,:,16,i) = V145(:,:,i)+dT.*dV145;
%     localinteractionbox(:,:,17,i) = V160(:,:,i)+dT.*dV160;
%     localinteractionbox(:,:,18,i) = V175(:,:,i)+dT.*dV175;
%     localinteractionbox(:,:,19,i) = V190(:,:,i)+dT.*dV190;
%     localinteractionbox(:,:,20,i) = V1105(:,:,i)+dT.*dV1105;
%     localinteractionbox(:,:,21,i) = V1120(:,:,i)+dT.*dV1120;
%     localinteractionbox(:,:,22,i) = V1135(:,:,i)+dT.*dV1135;
%     localinteractionbox(:,:,23,i) = V1150(:,:,i)+dT.*dV1150;
%     localinteractionbox(:,:,24,i) = V1165(:,:,i)+dT.*dV1165;
%     localinteractionbox(:,:,25,i) = V1(:,:,i)+dT.*dV1;
%     localinteractionbox(:,:,26,i) = V115(:,:,i)+dT.*dV115;
%     localinteractionbox(:,:,27,i) = V130(:,:,i)+dT.*dV130;
%     localinteractionbox(:,:,28,i) = V145(:,:,i)+dT.*dV145;
%     localinteractionbox(:,:,29,i) = V160(:,:,i)+dT.*dV160;
%     localinteractionbox(:,:,30,i) = V175(:,:,i)+dT.*dV175;
%     localinteractionbox(:,:,31,i) = V190(:,:,i)+dT.*dV190;
%     localinteractionbox(:,:,32,i) = V1105(:,:,i)+dT.*dV1105;
%     localinteractionbox(:,:,33,i) = V1120(:,:,i)+dT.*dV1120;
%     localinteractionbox(:,:,34,i) = V1135(:,:,i)+dT.*dV1135;
%     localinteractionbox(:,:,35,i) = V1150(:,:,i)+dT.*dV1150;
%     localinteractionbox(:,:,36,i) = V1165(:,:,i)+dT.*dV1165;
    
   localinteractionbox(:,:,1,i) = V190(:,:,i); % % four dimensions because three dimensional matrix + time
   
    localinteractionbox(:,:,2,i) = V1100(:,:,i);
   
    localinteractionbox(:,:,3,i) = V1110(:,:,i);
   
    localinteractionbox(:,:,4,i) = V1120(:,:,i);
    
    localinteractionbox(:,:,5,i) = V1130(:,:,i);
    
    localinteractionbox(:,:,6,i) = V1140(:,:,i);
   
    localinteractionbox(:,:,7,i) = V1150(:,:,i);
    
    localinteractionbox(:,:,8,i) = V1160(:,:,i);
   
    localinteractionbox(:,:,9,i) = V1170(:,:,i);
   
    localinteractionbox(:,:,10,i) = V1(:,:,i);
   
    localinteractionbox(:,:,11,i) = V110(:,:,i);
    
    localinteractionbox(:,:,12,i) = V120(:,:,i);
    
    localinteractionbox(:,:,13,i) = V130(:,:,i);
   
    localinteractionbox(:,:,14,i) = V140(:,:,i);
    
    localinteractionbox(:,:,15,i) = V150(:,:,i);
   
    localinteractionbox(:,:,16,i) = V160(:,:,i);
   
    localinteractionbox(:,:,17,i) = V170(:,:,i);
    
    localinteractionbox(:,:,18,i) = V180(:,:,i);
    
    localinteractionbox(:,:,19,i) = V190(:,:,i); % % four dimensions because three dimensional matrix + time
    
    localinteractionbox(:,:,20,i) = V1100(:,:,i);
  
    localinteractionbox(:,:,21,i) = V1110(:,:,i);
   
    localinteractionbox(:,:,22,i) = V1120(:,:,i);
   
    localinteractionbox(:,:,23,i) = V1130(:,:,i);
   
    localinteractionbox(:,:,24,i) = V1140(:,:,i);
    
    localinteractionbox(:,:,25,i) = V1150(:,:,i);
   
    localinteractionbox(:,:,26,i) = V1160(:,:,i);
   
    localinteractionbox(:,:,27,i) = V1170(:,:,i);
   
    localinteractionbox(:,:,28,i) = V1(:,:,i);
   
    localinteractionbox(:,:,29,i) = V110(:,:,i);
    
    localinteractionbox(:,:,30,i) = V120(:,:,i);
    
    localinteractionbox(:,:,31,i) = V130(:,:,i);
    
    localinteractionbox(:,:,32,i) = V140(:,:,i);
   
    localinteractionbox(:,:,33,i) = V150(:,:,i);
    
    localinteractionbox(:,:,34,i) = V160(:,:,i);
    
    localinteractionbox(:,:,35,i) = V170(:,:,i);
    
    localinteractionbox(:,:,36,i) = V180(:,:,i);
   
    
%              dlocal= -xa.*localinteraction(:,:,:,i)+...
%         1.*(1-localinteraction(:,:,:,i)).*is.*(convn(localinteractionbox(:,:,:,i), Interactione, 'same'))-...
%         1.*(1+localinteraction(:,:,:,i)).*is.*(convn(localinteractionbox(:,:,:,i), Interactioni, 'same'));
%         localinteraction(:,:,:,i+1)=localinteraction(:,:,:,i)+dT.*dlocal;

             dlocal= -xa.*localinteraction(:,:,:,i)+...
        1.*(1-localinteraction(:,:,:,i)).*is.*(convn(localinteractionbox(:,:,:,i), Interactione, 'same'))-...
        1.*(1+localinteraction(:,:,:,i)).*is.*(convn(localinteractionbox(:,:,:,i), Interactioni, 'same'));
        localinteraction(:,:,:,i+1)=localinteraction(:,:,:,i)+dT.*dlocal;
      
%         V2(:,:,i+1) = localinteraction(:,:,13,i)+dT.*dlocal(:,:,13);
%         V215(:,:,i+1) = localinteraction(:,:,14,i)+dT.*dlocal(:,:,14);
%         V230(:,:,i+1) = localinteraction(:,:,15,i)+dT.*dlocal(:,:,15);
%         V245(:,:,i+1) = localinteraction(:,:,16,i)+dT.*dlocal(:,:,16);
%         V260(:,:,i+1) = localinteraction(:,:,17,i)+dT.*dlocal(:,:,17);
%         V275(:,:,i+1) = localinteraction(:,:,18,i)+dT.*dlocal(:,:,18);
%         V290(:,:,i+1) = localinteraction(:,:,19,i)+dT.*dlocal(:,:,19);
%         V2105(:,:,i+1) = localinteraction(:,:,20,i)+dT.*dlocal(:,:,20);
%         V2120(:,:,i+1) = localinteraction(:,:,21,i)+dT.*dlocal(:,:,21);
%         V2135(:,:,i+1) = localinteraction(:,:,22,i)+dT.*dlocal(:,:,22);
%         V2150(:,:,i+1) = localinteraction(:,:,23,i)+dT.*dlocal(:,:,23);
%         V2165(:,:,i+1) = localinteraction(:,:,24,i)+dT.*dlocal(:,:,24);
        
%        V1feedback(:,:,i) = sum(localinteractionbox,3);

%      for w = 1:1:numwidth   %% segmentation of V1 neuronal field
%             endcol = w*squwidth;
%             begcol = endcol - (squwidth-1);
%             
%             read = V1(:,begcol:endcol,i);
%             read15 = V115(:,begcol:endcol,i);
%             read30 = V130(:,begcol:endcol,i);
%             read45 = V145(:,begcol:endcol,i);
%             read60 = V160(:,begcol:endcol,i);
%             read75 = V175(:,begcol:endcol,i);
%             read90 = V190(:,begcol:endcol,i);
%             read105 = V1105(:,begcol:endcol,i);
%             read120 = V1120(:,begcol:endcol,i);
%             read135 = V1135(:,begcol:endcol,i);
%             read150 =V1150(:,begcol:endcol,i);
%             read165 =V1165(:,begcol:endcol,i);
%             
% %             columnbox(:,1) = read.*sind(0);
% %             columnbox(:,2) = read15.*sind(15);
% %             columnbox(:,3) = read30.*sind(30);
% %             columnbox(:,4) = read45.*sind(45);
% %             columnbox(:,5) = read60.*sind(60);
% %             columnbox(:,6) = read75.*sind(75);
% %             columnbox(:,7) = read90.*sind(90);
% %             columnbox(:,8) = read105.*sind(105);
% %             columnbox(:,9) = read120.*sind(120);
% %             columnbox(:,10)=read135.*sind(135);
% %             columnbox(:,11)=read150.*sind(150);
% %             columnbox(:,12)=read165.*sind(165);
% 
%             columnbox(:,1) = read;
%             columnbox(:,2) = read15;
%             columnbox(:,3) = read30;
%             columnbox(:,4) = read45;
%             columnbox(:,5) = read60;
%             columnbox(:,6) = read75;
%             columnbox(:,7) = read90;
%             columnbox(:,8) = read105;
%             columnbox(:,9) = read120;
%             columnbox(:,10)=read135;
%             columnbox(:,11)=read150;
%             columnbox(:,12)=read165;
%             
%          dcol= -xa.*columnbox(:,:)+...
%         1.*(1-columnbox(:,:)).*(convn(columnbox(:,:), ee, 'same'))-...
%         1.*(1+columnbox(:,:)).*(convn(columnbox(:,:), iii, 'same'));
%         columnbox(:,:)=columnbox(:,:)+DT.*dcol;
%         
% % minVal = min(columnbox(:));
% % maxVal = max(columnbox(:));
% % columnbox = (columnbox - minVal) / ( maxVal - minVal );
% 
%         V1y(:,w) = columnbox(:,1);
%         V115y(:,w) = columnbox(:,2);
%         V130y(:,w) = columnbox(:,3);
%         V145y(:,w) = columnbox(:,4);
%         V160y(:,w) = columnbox(:,5);
%         V175y(:,w) = columnbox(:,6);
%         V190y(:,w) = columnbox(:,7);
%         V1105y(:,w) = columnbox(:,8);
%         V1120y(:,w) = columnbox(:,9);
%         V1135y(:,w) = columnbox(:,10);
%         V1150y(:,w) = columnbox(:,11);
%         V1165y(:,w) = columnbox(:,12);
%      end
%         for l = 1:1:numlength
%             % get the beginning and ending column for this square
%           
%             % get the beginning and ending row for this square
%             endrow = l*squlength;
%             begrow = endrow-(squlength-1);
%            
%             read = V1(begrow:endrow,:,i);
%             read15 = V115(begrow:endrow,:,i);
%             read30 = V130(begrow:endrow,:,i);
%             read45 = V145(begrow:endrow,:,i);
%             read60 = V160(begrow:endrow,:,i);
%             read75 = V175(begrow:endrow,:,i);
%             read90 = V190(begrow:endrow,:,i);
%             read105 = V1105(begrow:endrow,:,i);
%             read120 = V1120(begrow:endrow,:,i);
%             read135 = V1135(begrow:endrow,:,i);
%             read150 =V1150(begrow:endrow,:,i);
%             read165 =V1165(begrow:endrow,:,i);
%             
% %             rowbox(1,:) = read.*cosd(0);
% %             rowbox(2,:) = read15.*cosd(15);
% %             rowbox(3,:) = read30.*cosd(30);
% %             rowbox(4,:) = read45.*cosd(45);
% %             rowbox(5,:) = read60.*cosd(60);
% %             rowbox(6,:) = read75.*cosd(75);
% %             rowbox(7,:) = read90.*cosd(90);
% %             rowbox(8,:) = read105.*cosd(-15);
% %             rowbox(9,:) = read120.*cosd(-30);
% %             rowbox(10,:)=read135.*cosd(-45);
% %             rowbox(11,:)=read150.*cosd(-60);
% %             rowbox(12,:)=read165.*cosd(-75);
% 
%             rowbox(1,:) = read;
%             rowbox(2,:) = read15;
%             rowbox(3,:) = read30;
%             rowbox(4,:) = read45;
%             rowbox(5,:) = read60;
%             rowbox(6,:) = read75;
%             rowbox(7,:) = read90;
%             rowbox(8,:) = read105;
%             rowbox(9,:) = read120;
%             rowbox(10,:)=read135;
%             rowbox(11,:)=read150;
%             rowbox(12,:)=read165;
%           
%         drow= -xa.*rowbox(:,:)+...
%         1.*(1-rowbox(:,:)).*(convn(rowbox(:,:), ee, 'same'))-...
%         1.*(1+rowbox(:,:)).*(convn(rowbox(:,:), iii, 'same'));
%         rowbox(:,:)=rowbox(:,:)+DT.*drow;
% 
% % minVal1 = min(rowbox(:));
% % maxVal1 = max(rowbox(:));
% % rowbox = (rowbox - minVal1) / ( maxVal1 - minVal1 );
%         
%         V1x(l,:) = rowbox(1,:);
%         V115x(l,:) = rowbox(2,:);
%         V130x(l,:) = rowbox(3,:);
%         V145x(l,:) = rowbox(4,:);
%         V160x(l,:) = rowbox(5,:);
%         V175x(l,:) = rowbox(6,:);
%         V190x(l,:) = rowbox(7,:);
%         V1105x(l,:) = rowbox(8,:);
%         V1120x(l,:) = rowbox(9,:);
%         V1135x(l,:) = rowbox(10,:);
%         V1150x(l,:) = rowbox(11,:);
%         V1165x(l,:) = rowbox(12,:);
%         
%            negabox(negabox < 0) = -1;
%            negabox(negabox > 0) = 1;
%            negabox(negabox == 0) = 1;
%         end
%      
%      V1feedback(:,:,i) = V1x(:,:)+V1y(:,:)+V115x(:,:)+V115y(:,:)+V130x(:,:)+V130y(:,:)+V145x(:,:)+V145y(:,:)+V160x(:,:)+V160y(:,:)+V175x(:,:)+V175y(:,:)+V190x(:,:)+V190y(:,:)+V1105x(:,:)+V1105y(:,:)+V1120x(:,:)+V1120y(:,:)+V1135x(:,:)+V1135y(:,:)+V1150x(:,:)+V1150y(:,:)+V1165x(:,:)+V1165y(:,:);
% %      V1(:,:,i+1) = V1(:,:,i+1) + V1x(:,:)+V1y(:,:);
% %      V115(:,:,i+1) = V115(:,:,i+1) + V115x(:,:) + V1y(:,:);
%      V130(:,:,i+1) = V130(:,:,i+1) + V130x(:,:)+V130y(:,:);
%      V145(:,:,i+1) = V145(:,:,i+1) +V145x(:,:)+V145y(:,:);
%      V160(:,:,i+1) = V160(:,:,i+1) + V160x(:,:)+V160y(:,:);
%      V175(:,:,i+1)= V175(:,:,i+1)+V175x(:,:)+V175y(:,:);
%      V190(:,:,i+1)=V190(:,:,i+1)+V190x(:,:)+V190y(:,:);
%      V1105(:,:,i+1) = V1105(:,:,i+1)+V1105x(:,:)+V1105y(:,:);
%      V1120(:,:,i+1)=V1120(:,:,i+1)+V1120x(:,:)+V1120y(:,:);
%      V1135(:,:,i+1)=V1135(:,:,i+1)+V1135x(:,:)+V1135y(:,:);
%      V1150(:,:,i+1) = V1150(:,:,i+1)+V1150x(:,:)+V1150y(:,:);
%      V1165(:,:,i+1)=V1165(:,:,i+1)+V1165x(:,:)+V1165y(:,:);
%    
%     V1feedback(:,:,i) = V1(:,:,i) + V115(:,:,i) + V130(:,:,i) + V145(:,:,i) + V160(:,:,i) + V175(:,:,i)+V190(:,:,i)+V1105(:,:,i)+V1120(:,:,i)+V1135(:,:,i)+V1150(:,:,i)+V1165(:,:,i);
%     
%     Intraconnection = ((Horizontalscanbox.^2)+(Verticalscanbox.^2)).^0.5;
%     Intraconnection = Intraconnection.*negabox;
%     ICC = Intraconnection;
% %     ICC = Intraconnection./comparebox; % intracortical connectivity coeffcient
% %     
%      for w = 1:1:numwidth   %% intracortical feedback loop
%         for l = 1:1:numlength
%             % get the beginning and ending column for this square
%             endcol = w*squwidth;
%             begcol = endcol - (squwidth-1);
%             % get the beginning and ending row for this square
%             endrow = l*squlength;
%             begrow = endrow-(squlength-1);
%            
%             dear = (V1(begrow:endrow,begcol:endcol,i));
%             dear15 = (V115(begrow:endrow,begcol:endcol,i));
%             dear30 = (V130(begrow:endrow,begcol:endcol,i));
%             dear45 = (V145(begrow:endrow,begcol:endcol,i));
%             dear60 = (V160(begrow:endrow,begcol:endcol,i));
%             dear75 = (V175(begrow:endrow,begcol:endcol,i));
%             dear90 = (V190(begrow:endrow,begcol:endcol,i));
%             dear105 = (V1105(begrow:endrow,begcol:endcol,i));
%             dear120 = (V1120(begrow:endrow,begcol:endcol,i));
%             dear135 = (V1135(begrow:endrow,begcol:endcol,i));
%             dear150 =(V1150(begrow:endrow,begcol:endcol,i));
%             dear165 =(V1165(begrow:endrow,begcol:endcol,i));
%             
%             dear = dear.*ICC(l,w,1);
%             dear15=dear15.*ICC(l,w,2);
%             dear30=dear30.*ICC(l,w,3);
%             dear45=dear45.*ICC(l,w,4);
%             dear60=dear60.*ICC(l,w,5);
%             dear75=dear75.*ICC(l,w,6);
%             dear90=dear90.*ICC(l,w,7);
%             dear105=dear105.*ICC(l,w,8);
%             dear120=dear120.*ICC(l,w,9);
%             dear135=dear135.*ICC(l,w,10);
%             dear150=dear150.*ICC(l,w,11);
%             dear165=dear165.*ICC(l,w,12);
%             
%                  row1 = 1;
%                  column1 = 1;
% row2 = row1 + squlength - 1;
% column2 = column1 + squlength - 1;
% ro3 = squlength-1;
% co3 = squlength-1;
% 
% summatedlocalbox = dear+dear15+dear30+dear45+dear60+dear75+dear90+dear105+dear120+dear135+dear150+dear165;
% V1feedback((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = summatedlocalbox;  
%   V1((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear;  
%   V115((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear15;  
%   V130((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear30;  
%   V145((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear45; 
%   V160((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear60; 
%   V175((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear75; 
%   V190((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear90; 
%   V1105((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear105;
%   V1120((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear120; 
%   V1135((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear135; 
%   V1150((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear150; 
%   V1165((row2*l-ro3:row2*l), (column2*w-co3:column2*w),i) = dear165; 
%         end
%      end
%     V1feedback(:,:,i) = V1(:,:,i) + V115(:,:,i) + V130(:,:,i) + V145(:,:,i) + V160(:,:,i) + V175(:,:,i)+V190(:,:,i)+V1105(:,:,i)+V1120(:,:,i)+V1135(:,:,i)+V1150(:,:,i)+V1165(:,:,i);

%     V2090(:,:,i+1)=V1(:,:,i)+V190(:,:,i);
%     V215105(:,:,i+1)=V115(:,:,i)+V1105(:,:,i);
%     V230120(:,:,i+1)=V130(:,:,i)+V1120(:,:,i);
%     V245135(:,:,i+1)=V145(:,:,i)+V1135(:,:,i);
%     V260150(:,:,i+1)=V160(:,:,i)+V1150(:,:,i);
%     V275165(:,:,i+1)=V175(:,:,i)+V1165(:,:,i);
%     
    dV290= -xa.*V290(:,:,i)+...
        1.*(1-V290(:,:,i)).*fs.*(convn(localinteraction(:,:,19,i), V2verticale, 'same'))-...
        1.*(1+V290(:,:,i)).*fs.*(convn(localinteraction(:,:,19,i), V2verticali, 'same'));
    
    dV2= -xa.*V2(:,:,i)+...
        1.*(1-V2(:,:,i)).*fs.*(convn(localinteraction(:,:,10,i), V2horizontale, 'same'))-...
        1.*(1+V2(:,:,i)).*fs.*(convn(localinteraction(:,:,10,i), V2horizontali, 'same'));
    
    
        dV210= -xa.*V210(:,:,i)+...
        1.*(1-V210(:,:,i)).*fs.*(convn(localinteraction(:,:,11,i), V210e, 'same'))-...
        1.*(1+V210(:,:,i)).*fs.*(convn(localinteraction(:,:,11,i), V210i, 'same'));
    
            dV220= -xa.*V220(:,:,i)+...
        1.*(1-V220(:,:,i)).*fs.*(convn(localinteraction(:,:,12,i), V220e, 'same'))-...
        1.*(1+V220(:,:,i)).*fs.*(convn(localinteraction(:,:,12,i), V220i, 'same'));

    
            dV240= -xa.*V240(:,:,i)+...
        1.*(1-V240(:,:,i)).*fs.*(convn(localinteraction(:,:,14,i), V240e, 'same'))-...
        1.*(1+V240(:,:,i)).*fs.*(convn(localinteraction(:,:,14,i), V240i, 'same'));
    
    
            dV250= -xa.*V250(:,:,i)+...
        1.*(1-V250(:,:,i)).*fs.*(convn(localinteraction(:,:,15,i), V250e, 'same'))-...
        1.*(1+V250(:,:,i)).*fs.*(convn(localinteraction(:,:,15,i), V250i, 'same'));
    
    
    
            dV270= -xa.*V270(:,:,i)+...
        1.*(1-V270(:,:,i)).*fs.*(convn(localinteraction(:,:,17,i), V270e, 'same'))-...
        1.*(1+V270(:,:,i)).*fs.*(convn(localinteraction(:,:,17,i), V270i, 'same'));
    
            dV280= -xa.*V280(:,:,i)+...
        1.*(1-V280(:,:,i)).*fs.*(convn(localinteraction(:,:,18,i), V280e, 'same'))-...
        1.*(1+V280(:,:,i)).*fs.*(convn(localinteraction(:,:,18,i), V280i, 'same'));
    
    
            dV2100= -xa.*V2100(:,:,i)+...
        1.*(1-V2100(:,:,i)).*fs.*(convn(localinteraction(:,:,20,i), V2100e, 'same'))-...
        1.*(1+V2100(:,:,i)).*fs.*(convn(localinteraction(:,:,20,i), V2100i, 'same'));
    
            dV2110= -xa.*V2110(:,:,i)+...
        1.*(1-V2110(:,:,i)).*fs.*(convn(localinteraction(:,:,21,i), V2110e, 'same'))-...
        1.*(1+V2110(:,:,i)).*fs.*(convn(localinteraction(:,:,21,i), V2110i, 'same'));
    
    
            dV2130= -xa.*V2130(:,:,i)+...
        1.*(1-V2130(:,:,i)).*fs.*(convn(localinteraction(:,:,23,i), V2130e, 'same'))-...
        1.*(1+V2130(:,:,i)).*fs.*(convn(localinteraction(:,:,23,i), V2130i, 'same'));
    
            dV2140= -xa.*V2140(:,:,i)+...
        1.*(1-V2140(:,:,i)).*fs.*(convn(localinteraction(:,:,24,i), V2140e, 'same'))-...
        1.*(1+V2140(:,:,i)).*fs.*(convn(localinteraction(:,:,24,i), V2140i, 'same'));
    
    
            dV2160= -xa.*V2160(:,:,i)+...
        1.*(1-V2160(:,:,i)).*fs.*(convn(localinteraction(:,:,26,i), V2160e, 'same'))-...
        1.*(1+V2160(:,:,i)).*fs.*(convn(localinteraction(:,:,26,i), V2160i, 'same'));
    
            dV2170= -xa.*V2170(:,:,i)+...
        1.*(1-V2170(:,:,i)).*fs.*(convn(localinteraction(:,:,27,i), V2170e, 'same'))-...
        1.*(1+V2170(:,:,i)).*fs.*(convn(localinteraction(:,:,27,i), V2170i, 'same'));
    
    
    dV230= -xa.*V230(:,:,i)+...
        1.*(1-V230(:,:,i)).*fs.*(convn(localinteraction(:,:,13,i), V230e, 'same'))-...
        1.*(1+V230(:,:,i)).*fs.*(convn(localinteraction(:,:,13,i), V230i, 'same'));
    
    dV260= -xa.*V260(:,:,i)+...
        1.*(1-V260(:,:,i)).*fs.*(convn(localinteraction(:,:,16,i), V260e, 'same'))-...
        1.*(1+V260(:,:,i)).*fs.*(convn(localinteraction(:,:,16,i), V260i, 'same'));
    
    
    dV2120= -xa.*V2120(:,:,i)+...
        1.*(1-V2120(:,:,i)).*fs.*(convn(localinteraction(:,:,22,i), V2120e, 'same'))-...
        1.*(1+V2120(:,:,i)).*fs.*(convn(localinteraction(:,:,22,i), V2120i, 'same'));
    
    dV2150= -xa.*V2150(:,:,i)+...
        1.*(1-V2150(:,:,i)).*fs.*(convn(localinteraction(:,:,25,i), V2150e, 'same'))-...
        1.*(1+V2150(:,:,i)).*fs.*(convn(localinteraction(:,:,25,i), V2150i, 'same'));
    
    
    
    V2(:,:,i+1)=V2(:,:,i)+dT.*dV2;
  
    V210(:,:,i+1)=V210(:,:,i)+dT.*dV210;
   
    V220(:,:,i+1)=V220(:,:,i)+dT.*dV220;
   
    V230(:,:,i+1)=V230(:,:,i)+dT.*dV230;
  
    V240(:,:,i+1)=V240(:,:,i)+dT.*dV240;
   
    V250(:,:,i+1)=V250(:,:,i)+dT.*dV250;
 
    V260(:,:,i+1)=V260(:,:,i)+dT.*dV260;
   
    V270(:,:,i+1)=V270(:,:,i)+dT.*dV270;
 
    V280(:,:,i+1)=V280(:,:,i)+dT.*dV280;
  
    V290(:,:,i+1)=V290(:,:,i)+dT.*dV290;
   
    V2100(:,:,i+1)=V2100(:,:,i)+dT.*dV2100;
 
    V2110(:,:,i+1)=V2110(:,:,i)+dT.*dV2110;
  
    V2120(:,:,i+1)=V2120(:,:,i)+dT.*dV2120;
   
    V2130(:,:,i+1)=V2130(:,:,i)+dT.*dV2130;
   
    V2140(:,:,i+1)=V2140(:,:,i)+dT.*dV2140;
   
    V2150(:,:,i+1)=V2150(:,:,i)+dT.*dV2150;
   
    V2160(:,:,i+1)=V2160(:,:,i)+dT.*dV2160;
   
    V2170(:,:,i+1)=V2170(:,:,i)+dT.*dV2170;
    
    
% NormalizedV1(:,:,i+1)=V1+V122d5+V145+V167d5+V190+V1112d5+V1135+V1157d5;
end        % equation to model neuronal activity based on previous response
  
% compav = sum(Verticalscanbox,3); % %these are relics of my numerous tries
% compah = sum(Horizontalscanbox,3);
% compad = sum(meanbox1,3);
% exv = Verticalscanbox.^2;
% exh = Horizontalscanbox.^2;
% sums = (exv + exh).^0.5;
% df = sum(sums,3);
% sumc = meanbox1.^2+meanbox2.^2;
% sumc = sumc.^0.5;
% dc = sum(sumc,3);
% subplot(1,3,1); surf(compav); hold on; subplot(1,3,2); surf(compah); hold on; subplot(1,3,3); surf(compad)
% subplot(1,2,1); surf(df); hold on; subplot(1,2,2); surf(dc);
% tuning1 = [0,0.1,0.2,0.3,0.5,0.7,1,0.7,0.5,0.3,0.2,0.1];
% tuning2 = [0.3,0.5,0.7,1,0.7,0.5,0.3,0.2,0.1,0,0,0];
% angle=[0,15,30,45,60,75,90,105,120,135,150,165];

% xr = tuning1+tuning2+tuning1+tuning2;
% xs = xr.*0;
% xb = xr.*0;
% %%dT = 0.05; temporal delay affects the modulation rate of orientation
% %%tuning curve
% 
% % selective inhibition mechanism might have to do with vertical vs
% %%horizontal perception
% dV190= -xa.*xs(:,:)+...
%         1.*(1-xs(:,:)).*(convn(xr(:,:), e, 'same'))-...
%         1.*(1+xs(:,:)).*(convn(tuning2(:,:)+tuning2, ii, 'same'));
%     xs(:,:)=xs(:,:)+dT.*dV190;
    
%     vd = hypercolumn1 + hypercolumn2;
%     dV190= -xa.*xs(:,:)+...
%         1.*(1-xs(:,:)).*(convn(vd(:,:), e, 'same'))-...
%         1.*(1+xs(:,:)).*(convn(hypercolumn2(:,:), ii, 'same'));
%     xs(:,:)=xs(:,:)+dT.*dV190;
    
%     plot(angle,xs); hold on; plot(angle,tuning1); hold on; plot(angle,tuning2); xlabel('angle');legend('integrated tuning curve','local tuning curve 1','local tuning curve 2','fontsize',10,'location','northeast')
%     
%    dV190= -xa.*xb(:,:)+...
%         1.*(1-xb(:,:)).*(convn(xr(:,:), e, 'same'))-...
%         1.*(1+xb(:,:)).*(convn(tuning1(:,:), ii, 'same'));
%     xb(:,:)=xb(:,:)+dT.*dV190;
%      plot(angle,xb); hold on; plot(angle,tuning1); hold on; plot(angle,tuning2);xlabel('angle'); legend('integrated tuning curve','local tuning curve 1','local tuning curve 2','fontsize',10,'location','northeast')
%     

    % save output as a struct with mean and std as fields
%     output = struct;
%     output.mean = meanbox;
%     output.std = standard;
%     save(['grid' num2str(squwidth) '.mat'],'-struct','output')

figure(1);
surf(RGC(:,:,2));
% filename = 'smallreceptivefield.';
% saveas(figure(1),filename)

figure(2); % makes figure
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'animated.gif'; % creates file
for t = 1:size(V1feedback,3) % loops t over time steps. change 1st index to start at t != 0
    twodplotter = zeros(size(V1,1),size(V1,2)); % creates results array
    twodplotter(1:end,1:end) = V1feedback(:,:,t); % fills results array at t
    % contour3(twodplotter,size(RGC,3))
%     colormap(jet)
%     hold on
%     meshgrid(twodplotter,'Edgecolor', 'none');
    surf(twodplotter) % plots results
    title('V1') % names figure
    xlabel('Neuron Number (x)') % titles x axis, can change
    ylabel('Neuron Number (y)') % titles y axis, can change
    zlabel('Neural Activity') % titles z axis, can change
   % zlim(zlimit) % sets zlim based on zlimit
    drawnow 
      % Capture the plot as an image, which becomes one frame of the gif
      % Changes to figure appearance must occur before here
      frame = getframe(2);
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the gif file 
      if t == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end

figure(2); % makes figure
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'animated.gif'; % creates file
for t = 1:size(Interactione,3) % loops t over time steps. change 1st index to start at t != 0
    twodplotter = zeros(size(Interactione,1),size(Interactione,2)); % creates results array
    twodplotter(1:end,1:end) = Interactione(:,:,t); % fills results array at t
    % contour3(twodplotter,size(RGC,3))
%     colormap(jet)
%     hold on
%     meshgrid(twodplotter,'Edgecolor', 'none');
    surf(twodplotter) % plots results
    title('V1') % names figure
    xlabel('Neuron Number (x)') % titles x axis, can change
    ylabel('Neuron Number (y)') % titles y axis, can change
    zlabel('Neural Activity') % titles z axis, can change
   % zlim(zlimit) % sets zlim based on zlimit
    drawnow 
      % Capture the plot as an image, which becomes one frame of the gif
      % Changes to figure appearance must occur before here
      frame = getframe(2);
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the gif file 
      if t == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end

% Display the images created in subplots
hf = figure('units','normalized','position',[.2 .2 .6 .6]);
ax1 = subplot(2,3,1);
ibg = imagesc(Ie(:,:,1));
colormap(gca,lines)
axis off
title('Stimulus')
hcb1=colorbar;
ax2 = subplot(2,3,4);
iim = imagesc(V1(:,:,10));
colormap(gca,jet)
axis off
title('Neural representation')
hcb2=colorbar;
ax3 = subplot(2,3,[2:3, 5:6]);
ibg2 = surf(V1(:,:,10));
vmin = min(min(V1(:,:,10), [], 1:2));
vmax = max(max(V1(:,:,10), [], 1:2));
colormap(gca,jet);
caxis('manual')
% caxis(num2str(vmin), num2str(vmax));
% assignin('base','vmin',vmin); assignin('base','vmax',vmax)
caxis([vmin,vmax]);
axis off
hold on
iim2 = imagesc(Ie(:,:,1));
% colormap(gca,lines);
set(iim2,'AlphaData',0.3);
hcb3=colorbar;
title('Superimposed')


figure(4)
colormap(jet)
subplot(1,3,1); surf(RGC(:,:,10)); title('RGC');
hold on; subplot(1,3,2); surf(LGN(:,:,10)); title('LGN');
hold on; subplot(1,3,3); surf(V1(:,:,10)); title('V1');

figure(5)
colormap(jet)
subplot(1,3,1); imagesc(RGC(:,:,10)); title('RGC');
hold on; subplot(1,3,2); imagesc(LGN(:,:,10)); title('LGN');
hold on; subplot(1,3,3); imagesc(V1(:,:,10)); title('V1');

neuronrow = size(V2,1);
neuroncol = size(V2,2);
squwidth = 1;
squlength = squwidth;
verticalinaccurate = V2 +V210 + V220+ V230+V240 + +V2170+V2160++ V2150+V2140;
horizontalinaccurate = V290+V2100 + V2110+ V2120+V2130 + V250+ V260+V270+V280;
selectedx = V2;
selectedx(:,:,:) = 0;
selectedy = V2;
selectedy(:,:,:) = 0;

noverticalx(:,:,:) = (cosd(0).*V2(:,:,:))+(cosd(10).*V210(:,:,:))+(cosd(20).*V220(:,:,:))+(cosd(30).*V230(:,:,:))+(cosd(40).*V240(:,:,:))+(cosd(50).*V250(:,:,:))+(cosd(60).*V260(:,:,:))+(cosd(70).*V270(:,:,:))+(cosd(80).*V280(:,:,:))+(cosd(90).*V290(:,:,:))+(cosd(-80).*V2100(:,:,:))+(cosd(-70).*V2110(:,:,:))+(cosd(-60).*V2120(:,:,:))+(cosd(-50).*V2130(:,:,:))+(cosd(-40).*V2140(:,:,:))+(cosd(-30).*V2150(:,:,:))+(cosd(-20).*V2160(:,:,:))+(cosd(-10).*V2170(:,:,:));
nohorizontalx(:,:,:) = (cosd(0).*V2(:,:,:))+(cosd(10).*V210(:,:,:))+(cosd(20).*V220(:,:,:))+(cosd(30).*V230(:,:,:))+(cosd(40).*V240(:,:,:))+(cosd(50).*V250(:,:,:))+(cosd(60).*V260(:,:,:))+(cosd(70).*V270(:,:,:))+(cosd(80).*V280(:,:,:))+(cosd(90).*V290(:,:,:))+(cosd(100).*V2100(:,:,:))+(cosd(110).*V2110(:,:,:))+(cosd(120).*V2120(:,:,:))+(cosd(130).*V2130(:,:,:))+(cosd(140).*V2140(:,:,:))+(cosd(150).*V2150(:,:,:))+(cosd(160).*V2160(:,:,:))+(cosd(170).*V2170(:,:,:));
noverticaly(:,:,:) = (sind(0).*V2(:,:,:))+(sind(10).*V210(:,:,:))+(sind(20).*V220(:,:,:))+(sind(30).*V230(:,:,:))+(sind(40).*V240(:,:,:))+(sind(50).*V250(:,:,:))+(sind(60).*V260(:,:,:))+(sind(70).*V270(:,:,:))+(sind(80).*V280(:,:,:))+(sind(90).*V290(:,:,:))+(sind(-80).*V2100(:,:,:))+(sind(-70).*V2110(:,:,:))+(sind(-60).*V2120(:,:,:))+(sind(-50).*V2130(:,:,:))+(sind(-40).*V2140(:,:,:))+(sind(-30).*V2150(:,:,:))+(sind(-20).*V2160(:,:,:))+(sind(-10).*V2170(:,:,:));
nohorizontaly(:,:,:) = (sind(0).*V2(:,:,:))+(sind(10).*V210(:,:,:))+(sind(20).*V220(:,:,:))+(sind(30).*V230(:,:,:))+(sind(40).*V240(:,:,:))+(sind(50).*V250(:,:,:))+(sind(60).*V260(:,:,:))+(sind(70).*V270(:,:,:))+(sind(80).*V280(:,:,:))+(sind(90).*V290(:,:,:))+(sind(100).*V2100(:,:,:))+(sind(110).*V2110(:,:,:))+(sind(120).*V2120(:,:,:))+(sind(130).*V2130(:,:,:))+(sind(140).*V2140(:,:,:))+(sind(150).*V2150(:,:,:))+(sind(160).*V2160(:,:,:))+(sind(170).*V2170(:,:,:));

noverticalxx12 = zeros(size(V1,1),size(V1,2),size(noverticalx,3));
noverticalyy12 = zeros(size(V1,1),size(V1,2),size(noverticalx,3));
nohorizontalxx12 = zeros(size(V1,1),size(V1,2),size(noverticalx,3));
nohorizontalyy12= zeros(size(V1,1),size(V1,2),size(noverticalx,3));

for u = 1:size(vectorx)
noverticalxx12(vectorx(u),vectory(u),:) = noverticalx(vectorx(u),vectory(u),:);
noverticalyy12(vectorx(u),vectory(u),:) = noverticaly(vectorx(u),vectory(u),:);
nohorizontalxx12(vectorx(u),vectory(u),:) = nohorizontalx(vectorx(u),vectory(u),:);
nohorizontalyy12(vectorx(u),vectory(u),:) = nohorizontaly(vectorx(u),vectory(u),:);
end

for i = 1:1:neuroncol % this is the selection mechanism for orientation perception, for direction component, you need to also factor temporal component here
        for n = 1:1:neuronrow
            % get the beginning and ending column for this square
            endcol = i*squwidth;
            begcol = endcol - (squwidth-1);
            % get the beginning and ending row for this square
            endrow = n*squlength;
            begrow = endrow-(squlength-1);
            % get the data for that square
            read = abs(sum(verticalinaccurate(begrow:endrow,begcol:endcol,:)));
            read1 = abs(sum(horizontalinaccurate(begrow:endrow,begcol:endcol,:)));
            
            if read > read1
                selectedx(n,i,:) = noverticalxx12(n,i,:);
                selectedy(n,i,:) = noverticalyy12(n,i,:);
            elseif read < read1
                selectedx(n,i,:) = nohorizontalxx12(n,i,:);
                selectedy(n,i,:) = nohorizontalyy12(n,i,:);
            end

        end
end

s = 1; % arrow scaling factor
figure(3) % quivering
%subplot(1,2,1)
quivering00 = quiver(selectedx(:,:,time),selectedy(:,:,time),5);
quivering00.ShowArrowHead = 'off';
selectedvectorangles = atan2d(selectedx(:,:,time),selectedy(:,:,time));
hold on
iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% colormap(gca,lines);
set(iim234xxs,'AlphaData',0.2);
% hcb3334x=colorbar;
title('Selected Hypercolumn Response')

neuronrow = size(V1,1);
neuroncol = size(V1,2);
squwidth = 1;
squlength = squwidth;
verticalinaccurate1 = V1 +V110 +V120+ V130+V140 +V1170+V1160+ V1150+V1140;
horizontalinaccurate1 = V190+V1100 +V1110+ V1120+V1130 +V150+ V160+V170+V180;

selectedx1 = V1;
selectedx1(:,:,:) = 0;
selectedy1 = V1;
selectedy1(:,:,:) = 0;


noverticalx1(:,:,:) = (cosd(0).*V1(:,:,:))+(cosd(10).*V110(:,:,:))+(cosd(20).*V120(:,:,:))+(cosd(30).*V130(:,:,:))+(cosd(40).*V140(:,:,:))+(cosd(50).*V150(:,:,:))+(cosd(60).*V160(:,:,:))+(cosd(70).*V170(:,:,:))+(cosd(80).*V180(:,:,:))+(cosd(90).*V190(:,:,:))+(cosd(-80).*V1100(:,:,:))+(cosd(-70).*V1110(:,:,:))+(cosd(-60).*V1120(:,:,:))+(cosd(-50).*V1130(:,:,:))+(cosd(-40).*V1140(:,:,:))+(cosd(-30).*V1150(:,:,:))+(cosd(-20).*V1160(:,:,:))+(cosd(-10).*V1170(:,:,:));
nohorizontalx1(:,:,:) = (cosd(0).*V1(:,:,:))+(cosd(10).*V110(:,:,:))+(cosd(20).*V120(:,:,:))+(cosd(30).*V130(:,:,:))+(cosd(40).*V140(:,:,:))+(cosd(50).*V150(:,:,:))+(cosd(60).*V160(:,:,:))+(cosd(70).*V170(:,:,:))+(cosd(80).*V180(:,:,:))+(cosd(90).*V190(:,:,:))+(cosd(100).*V1100(:,:,:))+(cosd(110).*V1110(:,:,:))+(cosd(120).*V1120(:,:,:))+(cosd(130).*V1130(:,:,:))+(cosd(140).*V1140(:,:,:))+(cosd(150).*V1150(:,:,:))+(cosd(160).*V1160(:,:,:))+(cosd(170).*V1170(:,:,:));
noverticaly1(:,:,:) = (sind(0).*V1(:,:,:))+(sind(10).*V110(:,:,:))+(sind(20).*V120(:,:,:))+(sind(30).*V130(:,:,:))+(sind(40).*V140(:,:,:))+(sind(50).*V150(:,:,:))+(sind(60).*V160(:,:,:))+(sind(70).*V170(:,:,:))+(sind(80).*V180(:,:,:))+(sind(90).*V190(:,:,:))+(sind(-80).*V1100(:,:,:))+(sind(-70).*V1110(:,:,:))+(sind(-60).*V1120(:,:,:))+(sind(-50).*V1130(:,:,:))+(sind(-40).*V1140(:,:,:))+(sind(-30).*V1150(:,:,:))+(sind(-20).*V1160(:,:,:))+(sind(-10).*V1170(:,:,:));
nohorizontaly1(:,:,:) = (sind(0).*V1(:,:,:))+(sind(10).*V110(:,:,:))+(sind(20).*V120(:,:,:))+(sind(30).*V130(:,:,:))+(sind(40).*V140(:,:,:))+(sind(50).*V150(:,:,:))+(sind(60).*V160(:,:,:))+(sind(70).*V170(:,:,:))+(sind(80).*V180(:,:,:))+(sind(90).*V190(:,:,:))+(sind(100).*V1100(:,:,:))+(sind(110).*V1110(:,:,:))+(sind(120).*V1120(:,:,:))+(sind(130).*V1130(:,:,:))+(sind(140).*V1140(:,:,:))+(sind(150).*V1150(:,:,:))+(sind(160).*V1160(:,:,:))+(sind(170).*V1170(:,:,:));

noverticalxx1 = zeros(size(V1,1),size(V1,2),size(noverticalx1,3));
noverticalyy1 = zeros(size(V1,1),size(V1,2),size(noverticalx1,3));
nohorizontalxx1 = zeros(size(V1,1),size(V1,2),size(noverticalx1,3));
nohorizontalyy1= zeros(size(V1,1),size(V1,2),size(noverticalx1,3));

for u = 1:size(vectorx)
noverticalxx1(vectorx(u),vectory(u),:) = noverticalx1(vectorx(u),vectory(u),:);
noverticalyy1(vectorx(u),vectory(u),:) = noverticaly1(vectorx(u),vectory(u),:);
nohorizontalxx1(vectorx(u),vectory(u),:) = nohorizontalx1(vectorx(u),vectory(u),:);
nohorizontalyy1(vectorx(u),vectory(u),:) = nohorizontaly1(vectorx(u),vectory(u),:);
end

for i = 1:1:neuroncol % this is the selection mechanism for orientation perception, for direction component, you need to also factor temporal component here
        for n = 1:1:neuronrow
            % get the beginning and ending column for this square
            endcol = i*squwidth;
            begcol = endcol - (squwidth-1);
            % get the beginning and ending row for this square
            endrow = n*squlength;
            begrow = endrow-(squlength-1);
            % get the data for that square
            read = abs(sum(verticalinaccurate1(begrow:endrow,begcol:endcol,:)));
            read1 = abs(sum(horizontalinaccurate1(begrow:endrow,begcol:endcol,:)));
            
            if read > read1
                selectedx1(n,i,:) = noverticalxx1(n,i,:);
                selectedy1(n,i,:) = noverticalyy1(n,i,:);
            elseif read < read1
                selectedx1(n,i,:) = nohorizontalxx1(n,i,:);
                selectedy1(n,i,:) = nohorizontalyy1(n,i,:);
            end

        end
end

s = 1; % arrow scaling factor
 % quivering
%subplot(1,2,1)
quivering0 = quiver(selectedx1(:,:,time),selectedy1(:,:,time),5);
quivering0.ShowArrowHead = 'off';
selectedvectorangles11 = atan2d(selectedx1(:,:,time),selectedy1(:,:,time));
hold on
iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% colormap(gca,lines);
set(iim234xxs,'AlphaData',0.2);
% hcb3334x=colorbar;
title('Selected Hypercolumn Response')


    
neuronrow = size(V1,1);
neuroncol = size(V1,2);
squwidth = 1;
squlength = squwidth;
                                                                                                                                                                                                                         
                                                                                                                                                                                                                              
verticalinaccurate12 = localinteraction(:,:,10,:) +localinteraction(:,:,11,:)+localinteraction(:,:,12,:)+localinteraction(:,:,13,:)+ localinteraction(:,:,14,:)+localinteraction(:,:,27,:)+localinteraction(:,:,26,:)+localinteraction(:,:,25,:)+localinteraction(:,:,24,:);
horizontalinaccurate12 = localinteraction(:,:,19,:)+localinteraction(:,:,20,:)+localinteraction(:,:,21,:)+localinteraction(:,:,22,:)+localinteraction(:,:,23,:)+localinteraction(:,:,15,:)+localinteraction(:,:,16,:)+ localinteraction(:,:,17,:) + localinteraction(:,:,18,:);
selectedx12 = V1;
selectedx12(:,:,:) = 0;
selectedy12 = V1;
selectedy12(:,:,:) = 0;
noverticalxx = zeros(size(V1,1),size(V1,2),size(V1,3));
noverticalyy = zeros(size(V1,1),size(V1,2),size(V1,3));
nohorizontalxx = zeros(size(V1,1),size(V1,2),size(V1,3));
nohorizontalyy = zeros(size(V1,1),size(V1,2),size(V1,3));

 
noverticalx12(:,:,:) = (cosd(0).*localinteraction(:,:,10,:))+(cosd(10).*localinteraction(:,:,11,:))+(cosd(20).*localinteraction(:,:,12,:))+(cosd(30).*localinteraction(:,:,13,:))+(cosd(40).*localinteraction(:,:,14,:))+(cosd(50).*localinteraction(:,:,15,:))+(cosd(60).*localinteraction(:,:,16,:))+(cosd(70).*localinteraction(:,:,17,:))+(cosd(80).*localinteraction(:,:,18,:))+(cosd(90).*localinteraction(:,:,19,:))+(cosd(-80).*localinteraction(:,:,20,:))+(cosd(-70).*localinteraction(:,:,21,:))+(cosd(-60).*localinteraction(:,:,22,:))+(cosd(-50).*localinteraction(:,:,23,:))+(cosd(-40).*localinteraction(:,:,24,:))+(cosd(-30).*localinteraction(:,:,25,:))+(cosd(-20).*localinteraction(:,:,26,:))+(cosd(-10).*localinteraction(:,:,27,:));
nohorizontalx12(:,:,:) = (cosd(0).*localinteraction(:,:,10,:))+(cosd(10).*localinteraction(:,:,11,:))+(cosd(20).*localinteraction(:,:,12,:))+(cosd(30).*localinteraction(:,:,13,:))+(cosd(40).*localinteraction(:,:,14,:))+(cosd(50).*localinteraction(:,:,15,:))+(cosd(60).*localinteraction(:,:,16,:))+(cosd(70).*localinteraction(:,:,17,:))+(cosd(80).*localinteraction(:,:,18,:))+(cosd(90).*localinteraction(:,:,19,:))+(cosd(100).*localinteraction(:,:,20,:))+(cosd(110).*localinteraction(:,:,21,:))+(cosd(120).*localinteraction(:,:,22,:))+(cosd(130).*localinteraction(:,:,23,:))+(cosd(140).*localinteraction(:,:,24,:))+(cosd(150).*localinteraction(:,:,25,:))+(cosd(160).*localinteraction(:,:,26,:))+(cosd(170).*localinteraction(:,:,27,:));
noverticaly12(:,:,:) = (sind(0).*localinteraction(:,:,10,:))+(sind(10).*localinteraction(:,:,11,:))+(sind(20).*localinteraction(:,:,12,:))+(sind(30).*localinteraction(:,:,13,:))+(sind(40).*localinteraction(:,:,14,:))+(sind(50).*localinteraction(:,:,15,:))+(sind(60).*localinteraction(:,:,16,:))+(sind(70).*localinteraction(:,:,17,:))+(sind(80).*localinteraction(:,:,18,:))+(sind(90).*localinteraction(:,:,19,:))+(sind(-80).*localinteraction(:,:,20,:))+(sind(-70).*localinteraction(:,:,21,:))+(sind(-60).*localinteraction(:,:,22,:))+(sind(-50).*localinteraction(:,:,23,:))+(sind(-40).*localinteraction(:,:,24,:))+(sind(-30).*localinteraction(:,:,25,:))+(sind(-20).*localinteraction(:,:,26,:))+(sind(-10).*localinteraction(:,:,27,:));
nohorizontaly12(:,:,:) = (sind(0).*localinteraction(:,:,10,:))+(sind(10).*localinteraction(:,:,11,:))+(sind(20).*localinteraction(:,:,12,:))+(sind(30).*localinteraction(:,:,13,:))+(sind(40).*localinteraction(:,:,14,:))+(sind(50).*localinteraction(:,:,15,:))+(sind(60).*localinteraction(:,:,16,:))+(sind(70).*localinteraction(:,:,17,:))+(sind(80).*localinteraction(:,:,18,:))+(sind(90).*localinteraction(:,:,19,:))+(sind(100).*localinteraction(:,:,20,:))+(sind(110).*localinteraction(:,:,21,:))+(sind(120).*localinteraction(:,:,22,:))+(sind(130).*localinteraction(:,:,23,:))+(sind(140).*localinteraction(:,:,24,:))+(sind(150).*localinteraction(:,:,25,:))+(sind(160).*localinteraction(:,:,26,:))+(sind(170).*localinteraction(:,:,27,:));

for u = 1:size(vectorx)
noverticalxx(vectorx(u),vectory(u),:) = noverticalx12(vectorx(u),vectory(u),:);
noverticalyy(vectorx(u),vectory(u),:) = noverticaly12(vectorx(u),vectory(u),:);
nohorizontalxx(vectorx(u),vectory(u),:) = nohorizontalx12(vectorx(u),vectory(u),:);
nohorizontalyy(vectorx(u),vectory(u),:) = nohorizontaly12(vectorx(u),vectory(u),:);
end

for i = 1:1:neuroncol % this is the selection mechanism for orientation perception, for direction component, you need to also factor temporal component here
        for n = 1:1:neuronrow
            % get the beginning and ending column for this square
            endcol = i*squwidth;
            begcol = endcol - (squwidth-1);
            % get the beginning and ending row for this square
            endrow = n*squlength;
            begrow = endrow-(squlength-1);
            % get the data for that square
            read = abs(sum(verticalinaccurate12(begrow:endrow,begcol:endcol,:)));
            read1 = abs(sum(horizontalinaccurate12(begrow:endrow,begcol:endcol,:)));
            
            if read > read1
                selectedx12(n,i,:) = noverticalxx(n,i,:);
                selectedy12(n,i,:) = noverticalyy(n,i,:);
            elseif read < read1
                selectedx12(n,i,:) = nohorizontalxx(n,i,:);
                selectedy12(n,i,:) = nohorizontalyy(n,i,:);
            end

        end
end

s = 1; % arrow scaling factor
 % quivering
%subplot(1,2,1)
quivering0 = quiver(selectedx12(:,:,time),selectedy12(:,:,time),5);
quivering0.ShowArrowHead = 'off';
selectedvectorangles1 = atan2d(selectedx12(:,:,time),selectedy12(:,:,time));
hold on
iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% colormap(gca,lines);
set(iim234xxs,'AlphaData',0.2);
% hcb3334x=colorbar;
title('Selected Hypercolumn Response')


s = 1; % arrow scaling factor
figure(4) % quivering
%subplot(1,2,1)
quivering0 = quiver(noverticalx12(:,:,time),noverticaly12(:,:,time),5);
quivering0.ShowArrowHead = 'off';
noverticalvectorangles = atan2d(noverticalx(:,:,time),noverticaly(:,:,time));
hold on
iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% colormap(gca,lines);
set(iim234xxs,'AlphaData',0.2);
% hcb3334x=colorbar;
title('No Vertical Hypercolumn Response')

s = 1; % arrow scaling factor
figure(5) % quivering
%subplot(1,2,1)
quivering0 = quiver(nohorizontalx12(:,:,30),nohorizontaly12(:,:,30),5);
quivering0.ShowArrowHead = 'off';
nohorizontalvectorangles = atan2d(nohorizontalx(:,:,30),nohorizontaly(:,:,time));
hold on
iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% colormap(gca,lines);
set(iim234xxs,'AlphaData',0.2);
% hcb3334x=colorbar;
title('No Horizontal Hypercolumn Response')

subplot(2,2,1); df = quiver(noverticalx(:,:,30),noverticaly(:,:,30),5); df.ShowArrowHead = 'off'; title('novertical');
hold on; subplot(2,2,2); de = quiver(nohorizontalx(:,:,30),nohorizontaly(:,:,30),5); de.ShowArrowHead = 'off'; title('nohorizontal');
hold on; subplot(2,2,3); surf(noverticalvectorangles); 
hold on; subplot(2,2,4); surf(nohorizontalvectorangles)
% s = 1; % arrow scaling factor
% figure(6) % quivering
% %subplot(1,2,1)
% quivering0 = quiver((cosd(0).*V2(1:s:end,1:s:end,10)),(sind(0).*V2(1:s:end,1:s:end,10)),5);
% quivering0.ShowArrowHead = 'off';
% quiverx0 = quivering0.UData;
% quivery0 = quivering0.VData;
% vectorangle0 = atan2d(quivery0,quiverx0);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 0')
% 
% s = 1; % arrow scaling factor
% figure(7)
% quivering15 = quiver(cosd(15).*V215(1:s:end,1:s:end,10),sind(15).*V215(1:s:end,1:s:end,10),5);
% quivering15.ShowArrowHead = 'off';
% quiverx15 = quivering15.UData;
% quivery15 = quivering15.VData;
% vectorangle15 = atan2d(quivery15,quiverx15);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 15')
% 
% s = 1; % arrow scaling factor
% figure(8)
% quivering30 = quiver(cosd(30).*V230(1:s:end,1:s:end,10),sind(30).*V230(1:s:end,1:s:end,10),5);
% quivering30.ShowArrowHead = 'off';
% quiverx30 = quivering30.UData;
% quivery30 = quivering30.VData;
% vectorangle30 = atan2d(quivery30,quiverx30);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 30')
% 
% s = 1; % arrow scaling factor
% figure(9)
% quivering45 = quiver(cosd(45).*V245(1:s:end,1:s:end,10),sind(45).*V245(1:s:end,1:s:end,10),5);
% quivering45.ShowArrowHead = 'off';
% quiverx45 = quivering45.UData;
% quivery45 = quivering45.VData;
% vectorangle45 = atan2d(quivery45,quiverx45);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 45')
% 
% s = 1; % arrow scaling factor
% figure(10)
% quivering60 = quiver(cosd(60).*V260(1:s:end,1:s:end,10),sind(60).*V260(1:s:end,1:s:end,10),5);
% quivering60.ShowArrowHead = 'off';
% quiverx60 = quivering60.UData;
% quivery60 = quivering60.VData;
% vectorangle60 = atan2d(quivery60,quiverx60);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 60')
% 
% s = 1; % arrow scaling factor
% figure(11)
% quivering75 = quiver(cosd(75).*V275(1:s:end,1:s:end,10),sind(75).*V275(1:s:end,1:s:end,10),5);
% quivering75.ShowArrowHead = 'off';
% quiverx75 = quivering75.UData;
% quivery75 = quivering75.VData;
% vectorangle75 = atan2d(quivery75,quiverx75);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 75')
% 
% s = 1; % arrow scaling factor
% figure(12)
% quivering90 = quiver(cosd(90).*V290(1:s:end,1:s:end,10),sind(90).*V290(1:s:end,1:s:end,10),5);
% quivering90.ShowArrowHead = 'off';
% quiverx90 = quivering90.UData;
% quivery90 = quivering90.VData;
% vectorangle90 = atan2d(quivery90,quiverx90);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 90')
% 
% s = 1; % arrow scaling factor
% figure(13)
% quivering105 = quiver(cosd(105).*V2105(1:s:end,1:s:end,10),sind(105).*V2105(1:s:end,1:s:end,10),5);
% quivering105.ShowArrowHead = 'off';
% quiverx105 = quivering105.UData;
% quivery105 = quivering105.VData;
% vectorangle105 = atan2d(quivery105,quiverx105);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 105')
% 
% s = 1; % arrow scaling factor
% figure(14)
% quivering120 = quiver(cosd(120).*V2120(1:s:end,1:s:end,10),sind(120).*V2120(1:s:end,1:s:end,10),5);
% quivering120.ShowArrowHead = 'off';
% quiverx120 = quivering120.UData;
% quivery120 = quivering120.VData;
% vectorangle120 = atan2d(quivery120,quiverx120);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 120')
% 
% s = 1; % arrow scaling factor
% figure(15)
% quivering135 = quiver(cosd(135).*V2135(1:s:end,1:s:end,10),sind(135).*V2135(1:s:end,1:s:end,10),5);
% quivering135.ShowArrowHead = 'off';
% quiverx135 = quivering135.UData;
% quivery135 = quivering135.VData;
% vectorangle135 = atan2d(quivery135,quiverx135);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 135')
% 
% s = 1; % arrow scaling factor
% figure(16)
% quivering150 = quiver(cosd(150).*V2150(1:s:end,1:s:end,10),sind(150).*V2150(1:s:end,1:s:end,10),5);
% %quivering150.ShowArrowHead = 'off';
% quiverx150 = quivering150.UData;
% quivery150 = quivering150.VData;
% vectorangle150 = atan2d(quivery150,quiverx150);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 150')
% 
% s = 1; % arrow scaling factor
% figure(17)
% quivering165 = quiver(cosd(165).*V2165(1:s:end,1:s:end,10),sind(165).*V2165(1:s:end,1:s:end,10),5);
% quivering165.ShowArrowHead = 'off';
% quiverx165 = quivering165.UData;
% quivery165 = quivering165.VData;
% vectorangle165 = atan2d(quivery165,quiverx165);
% hold on
% iim234xxs = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxs,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering 165')
% 
% summatedvectorx = quiverx165+quiverx150+quiverx135+quiverx120+quiverx105+quiverx90+quiverx75+quiverx60+quiverx45+quiverx30+quiverx15+quiverx0;
% summatedvectory = quivery165+quivery150+quivery135+quivery120+quivery105+quivery90+quivery75+quivery60+quivery45+quivery30+quivery15+quivery0;
% 
% figure(18) % quivering normalized cell activity
% % subplot(1,2,1);
% quiveringgg = quiver(summatedvectorx,summatedvectory,5);
%   quiveringgg.ShowArrowHead = 'off';
%  quiverxns= quiveringgg.UData;
%  quiveryns = quiveringgg.VData;
%  vectorangles = atan2d(quiveryns,quiverxns);
% hold on
% iim234xxsg = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxsg,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering-summated vector response');
% 
% figure(19) % quivering normalized cell activity
% % subplot(1,2,1);
% quiveringgg = quiver(noverticalx(:,:,50),noverticaly(:,:,50),5);
%   quiveringgg.ShowArrowHead = 'off';
%  quiverxns= quiveringgg.UData;
%  quiveryns = quiveringgg.VData;
%  vectorangles = atan2d(quiveryns,quiverxns);
% hold on
% iim234xxsg = imagesc(imresize((Ie(:,:,1)),[floor(size(V1,1))/s floor(size(V1,2)/s)]));
% % colormap(gca,lines);
% set(iim234xxsg,'AlphaData',0.2);
% % hcb3334x=colorbar;
% title('Quivering-summated vector response');

% figure(22) % neuronal response from all orientations
% subplot(3,4,1)
% vvminv = min(min(V2(:,:,10), [], 1:2));
% vvmaxv = max(max(V2(:,:,10), [], 1:2));
% surf(V2(:,:,10)); colormap(gca,jet);
% caxis('manual')
% caxis([vvminv,vvmaxv]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% set(iim234xxxxv,'AlphaData',0.3);
% title('0 degree')
% subplot(3,4,2)
% surf(V215(:,:,10)); colormap(gca,jet);
% vvminvv = min(min(V215(:,:,10), [], 1:2));
% vvmaxvv = max(max(V215(:,:,10), [], 1:2));
% caxis('manual')
% caxis([vvminvv,vvmaxvv]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% set(iim234xxxxv,'AlphaData',0.3);
%     title('15 degree')
%     subplot(3,4,3)
% surf(V230(:,:,10)); colormap(jet);
% vvminvv30 = min(min(V230(:,:,10), [], 1:2));
% vvmaxvv30 = max(max(V230(:,:,10), [], 1:2));
% caxis('manual')
% caxis([vvminvv30,vvmaxvv30]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% % colormap(jet);
% set(iim234xxxxv,'AlphaData',0.3);
%     title('30 degree')
%        subplot(3,4,4)
% surf(V245(:,:,10)); colormap(jet);
% vvminvv45 = min(min(V245(:,:,10), [], 1:2));
% vvmaxvv45 = max(max(V245(:,:,10), [], 1:2));
% caxis('manual')
% caxis([vvminvv45,vvmaxvv45]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% set(iim234xxxxv,'AlphaData',0.3);
%     title('45 degree')
%       subplot(3,4,5)
% surf(V260(:,:,10)); colormap(jet);
% vvminvv60 = min(min(V260(:,:,10), [], 1:2));
% vvmaxvv60 = max(max(V260(:,:,10), [], 1:2));
% caxis('manual')
% caxis([vvminvv60,vvmaxvv60]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% set(iim234xxxxv,'AlphaData',0.3);
%     title('60 degree')
%      subplot(3,4,6)
% surf(V275(:,:,10)); colormap(jet);
% vvminvv75 = min(min(V275(:,:,10), [], 1:2));
% vvmaxvv75 = max(max(V275(:,:,10), [], 1:2));
% caxis('manual')
% caxis([vvminvv75,vvmaxvv75]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% set(iim234xxxxv,'AlphaData',0.3);
%     title('75 degree')
%        subplot(3,4,7)
% surf(V290(:,:,10)); colormap(jet);
% vvminvv90 = min(min(V290(:,:,10), [], 1:2));
% vvmaxvv90 = max(max(V290(:,:,10), [], 1:2));
% caxis('manual')
% caxis([vvminvv90,vvmaxvv90]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% set(iim234xxxxv,'AlphaData',0.3);
%     title('90 degree')
%      subplot(3,4,8)
% surf(V2105(:,:,10)); colormap(jet);
% vvminvv105 = min(min(V2105(:,:,10), [], 1:2));
% vvmaxvv105 = max(max(V2105(:,:,10), [], 1:2));
% caxis('manual')
% caxis([vvminvv105,vvmaxvv105]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% set(iim234xxxxv,'AlphaData',0.3);
%     title('105 degree')
%        subplot(3,4,9)
% surf(V2120(:,:,10)); colormap(jet);
% vvminvv120 = min(min(V2120(:,:,10), [], 1:2));
% vvmaxvv120 = max(max(V2120(:,:,10), [], 1:2));
% caxis('manual')
% caxis([vvminvv120,vvmaxvv120]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% set(iim234xxxxv,'AlphaData',0.3);
%     title('120 degree')
%        subplot(3,4,10)
% surf(V2135(:,:,10)); colormap(jet);
% vvminvv135 = min(min(V2135(:,:,10), [], 1:2));
% vvmaxvv135 = max(max(V2135(:,:,10), [], 1:2));
% caxis('manual')
% caxis([vvminvv135,vvmaxvv135]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% set(iim234xxxxv,'AlphaData',0.3);
%     title('135 degree')
%        subplot(3,4,11)
% surf(V2150(:,:,10)); colormap(jet);
% vvminvv150 = min(min(V2150(:,:,10), [], 1:2));
% vvmaxvv150 = max(max(V2150(:,:,10), [], 1:2));
% caxis('manual')
% caxis([vvminvv150,vvmaxvv150]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% set(iim234xxxxv,'AlphaData',0.3);
%     title('150 degree')
%        subplot(3,4,12)
% surf(V2165(:,:,10)); colormap(jet);
% vvminvv165 = min(min(V2165(:,:,10), [], 1:2));
% vvmaxvv165 = max(max(V2165(:,:,10), [], 1:2));
% caxis('manual')
% caxis([vvminvv165,vvmaxvv165]);
% axis off
% hold on
% iim234xxxxv = imagesc(Ie(:,:,1));
% set(iim234xxxxv,'AlphaData',0.3);
%     title('165 degree')
    


% a0=max(V1(:,:,10),[],1:2);
% a15=max(V115(:,:,10),[],1:2);
% a30=max(V130(:,:,10),[],1:2);
% a45=max(V145(:,:,10),[],1:2);
% a60=max(V160(:,:,10),[],1:2);
% a75=max(V175(:,:,10),[],1:2);
% a90=max(V190(:,:,10),[],1:2);
% a105=max(V1105(:,:,10),[],1:2);
% a120=max(V1120(:,:,10),[],1:2);
% a135=max(V1135(:,:,10),[],1:2);
% a150=max(V1150(:,:,10),[],1:2);
% a165=max(V1165(:,:,10),[],1:2);

% a0=mean(V1(:,:,10),1:2);
% a15=mean(V115(:,:,10),1:2);
% a30=mean(V130(:,:,10),1:2);
% a45=mean(V145(:,:,10),1:2);
% a60=mean(V160(:,:,10),1:2);
% a75=mean(V175(:,:,10),1:2);
% a90=mean(V190(:,:,10),1:2);
% a105=mean(V1105(:,:,10),1:2);
% a120=mean(V1120(:,:,10),1:2);
% a135=mean(V1135(:,:,10),1:2);
% a150=mean(V1150(:,:,10),1:2);
% a165=mean(V1165(:,:,10),1:2);

% a0=max(V1(:,:,10),[],1:2);
% a15=max(V115(:,:,10),[],1:2);
% a30=max(V130(:,:,10),[],1:2);
% a45=max(V145(:,:,10),[],1:2);
% a60=max(V160(:,:,10),[],1:2);
% a75=max(V175(:,:,10),[],1:2);
% a90=max(V190(:,:,10),[],1:2);
% a105=max(V1105(:,:,10),[],1:2);
% a120=max(V1120(:,:,10),[],1:2);
% a135=max(V1135(:,:,10),[],1:2);
% a150=max(V1150(:,:,10),[],1:2);
% a165=max(V1165(:,:,10),[],1:2);

a0=max(V2(:,:,50),[],1:2);
a15=max(V215(:,:,50),[],1:2);
a30=max(V230(:,:,50),[],1:2);
a45=max(V245(:,:,50),[],1:2);
a60=max(V260(:,:,50),[],1:2);
a75=max(V275(:,:,50),[],1:2);
a90=max(V290(:,:,50),[],1:2);
a105=max(V2105(:,:,50),[],1:2);
a120=max(V2120(:,:,50),[],1:2);
a135=max(V2135(:,:,50),[],1:2);
a150=max(V2150(:,:,50),[],1:2);
a165=max(V2165(:,:,50),[],1:2);

b0=max(V1(:,:,50),[],1:2);
b15=max(V115(:,:,50),[],1:2);
b30=max(V130(:,:,50),[],1:2);
b45=max(V145(:,:,50),[],1:2);
b60=max(V160(:,:,50),[],1:2);
b75=max(V175(:,:,50),[],1:2);
b90=max(V190(:,:,50),[],1:2);
b105=max(V1105(:,:,50),[],1:2);
b120=max(V1120(:,:,50),[],1:2);
b135=max(V1135(:,:,50),[],1:2);
b150=max(V1150(:,:,50),[],1:2);
b165=max(V1165(:,:,50),[],1:2);

c0 = max(V2(112:146,125:145,50),[],1:2);
c15 = max(V215(112:146,125:145,50),[],1:2);
c30 = max(V230(112:146,125:145,50),[],1:2);
c45 = max(V245(112:146,125:145,50),[],1:2);
c60 = max(V260(112:146,125:145,50),[],1:2);
c75 = max(V275(112:146,125:145,50),[],1:2);
c90 = max(V290(112:146,125:145,50),[],1:2);
c105 = max(V2105(112:146,125:145,50),[],1:2);
c120 = max(V2120(112:146,125:145,50),[],1:2);
c135 = max(V2135(112:146,125:145,50),[],1:2);
c150 = max(V2150(112:146,125:145,50),[],1:2);
c165 = max(V2165(112:146,125:145,50),[],1:2);

d0 = max(V2(125,135,50),[],1:2);
d15 = max(V215(125,135,50),[],1:2);
d30 = max(V230(125,135,50),[],1:2);
d45 = max(V245(125,135,50),[],1:2);
d60 = max(V260(125,135,50),[],1:2);
d75 = max(V275(125,135,50),[],1:2);
d90 = max(V290(125,135,50),[],1:2);
d105 = max(V2105(125,135,50),[],1:2);
d120 = max(V2120(125,135,50),[],1:2);
d135 = max(V2135(125,135,50),[],1:2);
d150 = max(V2150(125,135,50),[],1:2);
d165 = max(V2165(125,135,50),[],1:2);

% py = [a105,a120,a135,a150,a165,a0,a15,a30,a45,a60,a75,a90];

py = [a0,a15,a30,a45,a60,a75,a90,a105,a120,a135,a150,a165,a0];
by = [b0,b15,b30,b45,b60,b75,b90,b105,b120,b135,b150,b165,b0];
cy = [c0,c15,c30,c45,c60,c75,c90,c105,c120,c135,c150,c165,c0];
dy = [d0,d15,d30,d45,d60,d75,d90,d105,d120,d135,d150,d165,d0];
px = [0,15,30,45,60,75,90,105,120,135,150,165,180];

plot(px,py);
xlabel('Angle for which orientation cells are tuned','Fontsize',16);
ylabel('Neural Activation Level of Orientation Cells','Fontsize',16);
title('Orientation Response Competition','Fontsize',16);

plot(px,by);
xlabel('Angle for which orientation cells are tuned','Fontsize',16);
ylabel('Neural Activation Level of Orientation Cells','Fontsize',16);
title('Orientation Response Competition','Fontsize',16);

plot(px,cy);
xlabel('Angle for which orientation cells are tuned','Fontsize',16);
ylabel('Neural Activation Level of Orientation Cells','Fontsize',16);
title('Orientation Response Competition','Fontsize',16);

plot(px,dy); %% %% response to the vertical line only
xlabel('Angle for which orientation cells are tuned','Fontsize',16);
ylabel('Neural Activation Level of Orientation Cells','Fontsize',16);
title('Orientation Response Competition','Fontsize',16);

subplot(1,2,1);
plot(px,by);
xlabel('Angle for which orientation cells are tuned','Fontsize',11);
ylabel('Neural Activation Level of Orientation Cells','Fontsize',11);
title('before interactionfield','Fontsize',12);
hold on;
subplot(1,2,2);
plot(px,py);
xlabel('Angle for which orientation cells are tuned','Fontsize',11);
ylabel('Neural Activation Level of Orientation Cells','Fontsize',11);
title('after interactionfield','Fontsize',12);

%Interactionvisual = Interactione - Interactioni;

% for f = 1:size(Interactione,3)
%     Interactionvisual(:,:,f) = Interactione(:,:,f)*Interactioni(:,:,f);
%   %  Interactionvisual(:,f,:) = squeeze(Interactione(:,f,:))*squeeze(Interactioni(:,f,:));
% end
% c = fspecial3('gaussian',[40 40 40],[3 3 1]);
% d = fspecial3('gaussian',[40 40 40],[6 6 2]);
[X,Y,Z] = meshgrid(1:size(Interactioni,1), 1:size(Interactioni,2), 1:size(Interactioni,3));

Xc = X(:);
Yc = Y(:);
Zc = Z(:);
Datac = Interactioni(:);
downsample = 1; % displaying 5 million points makes MATLAB unhappy
scatter3(Xc(1:downsample:end),Yc(1:downsample:end),Zc(1:downsample:end),15,Interactioni(1:downsample:end),'filled')
colormap(lines)
vvminvv165 = (min(Interactioni(:,:,:), [], 1:3));
vvmaxvv165 = (max(Interactioni(:,:,:), [], 1:3));
caxis('manual')
caxis([vvminvv165,vvmaxvv165]);% default is usually jet(64)
colorbar


% [X,Y,Z] = sphere;
% surf(X,Y,Z);
% 
% 
% sigma = 1.76;
% %Window size
% sz = 4;
% [x,y,z]=meshgrid(-sz:sz,-sz:sz,-sz:sz);
% 
% M = size(x,1)-1;
% N = size(y,1)-1;
% P = size(z,1)-1;
% Exp_comp = -(x.^2+y.^2+z.^2)/(2*sigma*sigma);
% Kernel= exp(Exp_comp)/(2*pi*sigma*sigma);
% [x,y,z]=sphere;
% surf(x,y,z)


 hypercolumnslocal = localinteraction(:,:,10:28,time); % % individual hypercolumn responses after local interaction
 hypercolumnboxl = zeros(size(sortx,1),19);
hypercolumnl = zeros(1,19);
roots = ceil(size(sortx,1)^0.5);
half = ceil(size(sortx,1)*0.5); % % for poggendorff illusion
root = ceil(half^0.5);
subplotx = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180];

for s= 2:half
    hypercolumnl(:,:)=hypercolumnslocal(sortx(s),sorty(s),:);
    disp(hypercolumnl)
            % let's plot it
%            subplot(numlength,numwidth,(n-1)*numwidth+i)
%            plot(booksqueeze)
%            title([num2str(i) ' ' num2str(n)])
            hypercolumnboxl(s,:)=hypercolumnl;
             subplot(root,root,s)
             plot(subplotx,hypercolumnl); title(['x' num2str(sortx(s)) ' '  'y' num2str(sorty(s))]);
             hold on; 
             plot([1 1]*(90), ylim, '--k')
             hold off;
             output = struct;
    output.hypercolumns = hypercolumnboxl;
            save(['individualhypercolumns' 'interactionfieldsize' num2str(xsigmae) 'eiratio' num2str(eiratio) '.mat'],'-struct','output')
end

 hypercolumnslocal = localinteraction(:,:,10:28,time); % % individual hypercolumn responses after local interaction
 hypercolumnbox1l = zeros(size(vectorx,1),19);
hypercolumnl1 = zeros(1,19);
roots = ceil(size(vectorx,1)^0.5);
half = ceil(size(vectorx,1)*0.5); % % for poggendorff illusion
root = ceil(half^0.5);
subplotx = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180];

for s= 2:half
    hypercolumnl1(:,:)=hypercolumnslocal(vectorx(s),vectory(s),:);
    disp(hypercolumnl1)
            % let's plot it
%            subplot(numlength,numwidth,(n-1)*numwidth+i)
%            plot(booksqueeze)
%            title([num2str(i) ' ' num2str(n)])
            hypercolumnbox1l(s,:)=hypercolumn11;
             subplot(root,root,s)
             plot(subplotx,hypercolumnl1); title(['x' num2str(vectorx(s)) ' '  'y' num2str(vectory(s))]);
             hold on; 
             plot([1 1]*(110), ylim, '--k')
             hold off;
             
end



plot(subplotx,d); xlabel('Orientation (deg)','fontsize',16); ylabel('Neuronal Activation Level','fontsize',16);
hold on; 
plot([1 1]*(70), ylim, '--k')
hold off;

plot(subplotx,c);xlabel('Orientation (deg)','fontsize',16); ylabel('Neuronal Activation Level','fontsize',16);
hold on; 
plot([1 1]*(110), ylim, '--k')
hold off;


d = [-1.76462299405978e-07;-2.99177916105349e-07;-3.51605260308400e-07;-1.44244474105000e-07;4.77445873835882e-07;1.59601042549726e-06;3.14237183411158e-06;4.74606399535912e-06;5.74007036549980e-06;5.50503129011053e-06;3.98302359360677e-06;1.83036600209817e-06;-3.00919855125682e-08;-1.02852408566401e-06;-1.17769006031531e-06;-8.46410575919837e-07;-4.40527628546284e-07;-1.90537188641983e-07;-1.42386799286793e-07];
c = [-2.57889582478009e-07;-3.26610708528148e-07;-3.40298620942311e-07;-3.29969681270090e-07;-2.67804511412084e-07;-6.64335361335916e-08;4.23048958405924e-07;1.38817535790208e-06;2.91808772161046e-06;4.79380865812147e-06;6.41229294455187e-06;7.05897014946583e-06;6.40994628878628e-06;4.79212936790140e-06;2.92005006160530e-06;1.39711091490919e-06;4.40472898350020e-07;-3.59205586370883e-08;-2.13102844079085e-07];   

plot(subplotx,e); xlabel('Orientation (deg)','fontsize',16); ylabel('Neuronal Activation Level','fontsize',16);
hold on; 
plot([1 1]*(90), ylim, '--k')
hold off;

plot(a); hold on; plot(b); hold on; plot(c); xlim([0 121]); title('Receptive Field','fontsize',16); legend(' A. i/e ratio 1.1',' B. i/e ratio 1.6',' C. i/e ratio 2.1');
 
f = normalize(d);
g = normalize(c);
h = normalize(e);
% 
 plot(subplotx,f); hold on; plot(subplotx,g); hold on; plot(subplotx,h); xlabel('Orientation (deg)','fontsize',16); ylabel('Normalized Scale','fontsize',16); legend('Distorted perception','Tilted line','Vertical line','fontsize',16);

hypercolumns(:,:,1) = V2(:,:,time);% % individual hypercolumn responses of V2 after local interaction 
hypercolumns(:,:,2) = V210(:,:,time);
hypercolumns(:,:,3) = V220(:,:,time);
hypercolumns(:,:,4) = V230(:,:,time);
hypercolumns(:,:,5) = V240(:,:,time);
hypercolumns(:,:,6) = V250(:,:,time);
hypercolumns(:,:,7) = V260(:,:,time);
hypercolumns(:,:,8) = V270(:,:,time);
hypercolumns(:,:,9) = V280(:,:,time);
hypercolumns(:,:,10) = V290(:,:,time);
hypercolumns(:,:,11) = V2100(:,:,time);
hypercolumns(:,:,12) = V2110(:,:,time);
hypercolumns(:,:,13) = V2120(:,:,time);
hypercolumns(:,:,14) = V2130(:,:,time);
hypercolumns(:,:,15) = V2140(:,:,time);
hypercolumns(:,:,16) = V2150(:,:,time);
hypercolumns(:,:,17) = V2160(:,:,time);
hypercolumns(:,:,18) = V2170(:,:,time);
hypercolumns(:,:,19) = V2(:,:,time);
hypercolumnbox = zeros(size(sortx,1),19);
hypercolumn = zeros(1,19);
sidehypercolumns = localinteraction(:,:,10:18,time);
sidehypercolumnbox = zeros(size(sortx,1),19);
sidehypercolumn = zeros(1,19);
roots = ceil(size(sortx,1)^0.5);
half = ceil(size(sortx,1)*0.5); % % for poggendorff illusion
root = ceil(half^0.5);
subplotx = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180];
for s= 2:half
    hypercolumn(:,:)=hypercolumns(sortx(s),sorty(s),:);
    disp(hypercolumn)
            % let's plot it
%            subplot(numlength,numwidth,(n-1)*numwidth+i)
%            plot(booksqueeze)
%            title([num2str(i) ' ' num2str(n)])
            hypercolumnbox(s,:)=hypercolumn;
             subplot(root,root,s)
             plot(subplotx,hypercolumn); title(['x' num2str(sortx(s)) ' '  'y' num2str(sorty(s))]);
             hold on; 
             plot([1 1]*(90-localtheta), ylim, '--k')
             hold off;
             
end

for s=1:size(sidesortx,1)
               sidehypercolumn(:,:)=sidehypercolumns(sidesortx(s),sidesorty(s),:); % individual hypercolumn responses after local interaction
    disp(sidehypercolumn)
            % let's plot it
%            subplot(numlength,numwidth,(n-1)*numwidth+i)
%            plot(booksqueeze)
%            title([num2str(i) ' ' num2str(n)])
            sidehypercolumnbox(s,:)=sidehypercolumn;
%             standardbox = double(read);
%             standard(n,i,:) = std(standardbox,0,[1 2]);
end

hypercolumnsv1(:,:,1) = V1(:,:,time);% % individual hypercolumn responses of V2 after local interaction 
hypercolumnsv1(:,:,2) = V110(:,:,time);
hypercolumnsv1(:,:,3) = V120(:,:,time);
hypercolumnsv1(:,:,4) = V130(:,:,time);
hypercolumnsv1(:,:,5) = V140(:,:,time);
hypercolumnsv1(:,:,6) = V150(:,:,time);
hypercolumnsv1(:,:,7) = V160(:,:,time);
hypercolumnsv1(:,:,8) = V170(:,:,time);
hypercolumnsv1(:,:,9) = V180(:,:,time);
hypercolumnsv1(:,:,10) = V190(:,:,time);
hypercolumnsv1(:,:,11) = V1100(:,:,time);
hypercolumnsv1(:,:,12) = V1110(:,:,time);
hypercolumnsv1(:,:,13) = V1120(:,:,time);
hypercolumnsv1(:,:,14) = V1130(:,:,time);
hypercolumnsv1(:,:,15) = V1140(:,:,time);
hypercolumnsv1(:,:,16) = V1150(:,:,time);
hypercolumnsv1(:,:,17) = V1160(:,:,time);
hypercolumnsv1(:,:,18) = V1170(:,:,time);
hypercolumnsv1(:,:,19) = V1(:,:,time);
hypercolumnboxv1 = zeros(size(sortx,1),19);
hypercolumnv1 = zeros(1,19);
% sidehypercolumnsv1 = localinteractionv1(:,:,13:24,50);
% sidehypercolumnboxv1 = zeros(size(sortx,1),12);
% sidehypercolumnv1 = zeros(1,12);
roots = ceil(size(sortx,1)^0.5);
half = ceil(size(sortx,1)*0.5); % % for poggendorff illusion
root = ceil(half^0.5);
subplotx = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180];
for s= 2:half
    hypercolumnv1(:,:)=hypercolumnsv1(sortx(s),sorty(s),:);
    disp(hypercolumnv1)
            % let's plot it
%            subplot(numlength,numwidth,(n-1)*numwidth+i)
%            plot(booksqueeze)
%            title([num2str(i) ' ' num2str(n)])
            hypercolumnboxv1(s,:)=hypercolumnv1;
             subplot(root,root,s); 
             plot(subplotx,hypercolumnv1); title(['x' num2str(sortx(s)) ' '  'y' num2str(sorty(s))])
             hold on; 
             plot([1 1]*(90-localtheta), ylim, '--k')
             hold off;
             
             
end

% % % let's try investigating sub additivity of normalization
 subadd = V2+V215+V230+V245+V260+V275+V290+V2105+V2120+V2135+V2150+V2165;
  subadd = V1+V115+V130+V145+V160+V175+V190+V1105+V1120+V1135+V1150+V1165;
subadditivity = sum(sum(subadd,3),1:2);


[X,Y,Z] = meshgrid(1:size(Interactioni,1), 1:size(Interactioni,2), 1:size(Interactioni,3));

Xc = X(:);
Yc = Y(:);
Zc = Z(:);
Datac = Interactioni(:);
downsample = 1; % displaying 5 million points makes MATLAB unhappy
scatter3(Xc(1:downsample:end),Yc(1:downsample:end),Zc(1:downsample:end),15,Interactioni(1:downsample:end),'filled')
%colormap(jet
colormap(flipud(gray(4)));
vvminvv165 = (min(Interactioni(:,:,:), [], 1:3));
vvmaxvv165 = (max(Interactioni(:,:,:), [], 1:3));
caxis('manual')
caxis([vvminvv165,vvmaxvv165]);% default is usually jet(64)
% cmap = jet(size(Interactioni,1)) ;
% %%Arrange the colors range
% colors = zeros(size(Z));         
% for i = 1:size(Interactioni,1)
%     colors(Z > Interactioni(i,1) & Z<=Interactioni(i,2)) = Interactioni(i,2);           
% end
%colorbar

neuronrow = size(V1,1);
neuroncol = size(V1,2);
squwidth = 1;
squlength = squwidth;
                                                                                                                                                                                                                         
                                                                                                                                                                                                                              
verticalinaccurate12t = localinteraction(:,:,10,:) +localinteraction(:,:,11,:)+localinteraction(:,:,12,:)+localinteraction(:,:,13,:)+ localinteraction(:,:,14,:)+localinteraction(:,:,27,:)+localinteraction(:,:,26,:)+localinteraction(:,:,25,:)+localinteraction(:,:,24,:);
horizontalinaccurate12t = localinteraction(:,:,19,:)+localinteraction(:,:,20,:)+localinteraction(:,:,21,:)+localinteraction(:,:,22,:)+localinteraction(:,:,23,:)+localinteraction(:,:,15,:)+localinteraction(:,:,16,:)+ localinteraction(:,:,17,:) + localinteraction(:,:,18,:);
selectedx12t = V1;
selectedx12t(:,:,:) = 0;
selectedy12t = V1;
selectedy12t(:,:,:) = 0;
noverticalxxt = zeros(size(V1,1),size(V1,2),size(V1,3));
noverticalyyt = zeros(size(V1,1),size(V1,2),size(V1,3));
nohorizontalxxt = zeros(size(V1,1),size(V1,2),size(V1,3));
nohorizontalyyt = zeros(size(V1,1),size(V1,2),size(V1,3));
vectorangles = zeros(size(V1,1),size(V1,2),size(V1,3));


noverticalx12t(:,:,:) = (cosd(0).*localinteraction(:,:,10,:))+(cosd(10).*localinteraction(:,:,11,:))+(cosd(20).*localinteraction(:,:,12,:))+(cosd(30).*localinteraction(:,:,13,:))+(cosd(40).*localinteraction(:,:,14,:))+(cosd(50).*localinteraction(:,:,15,:))+(cosd(60).*localinteraction(:,:,16,:))+(cosd(70).*localinteraction(:,:,17,:))+(cosd(80).*localinteraction(:,:,18,:))+(cosd(90).*localinteraction(:,:,19,:))+(cosd(-80).*localinteraction(:,:,20,:))+(cosd(-70).*localinteraction(:,:,21,:))+(cosd(-60).*localinteraction(:,:,22,:))+(cosd(-50).*localinteraction(:,:,23,:))+(cosd(-40).*localinteraction(:,:,24,:))+(cosd(-30).*localinteraction(:,:,25,:))+(cosd(-20).*localinteraction(:,:,26,:))+(cosd(-10).*localinteraction(:,:,27,:));
nohorizontalx12t(:,:,:) = (cosd(0).*localinteraction(:,:,10,:))+(cosd(10).*localinteraction(:,:,11,:))+(cosd(20).*localinteraction(:,:,12,:))+(cosd(30).*localinteraction(:,:,13,:))+(cosd(40).*localinteraction(:,:,14,:))+(cosd(50).*localinteraction(:,:,15,:))+(cosd(60).*localinteraction(:,:,16,:))+(cosd(70).*localinteraction(:,:,17,:))+(cosd(80).*localinteraction(:,:,18,:))+(cosd(90).*localinteraction(:,:,19,:))+(cosd(100).*localinteraction(:,:,20,:))+(cosd(110).*localinteraction(:,:,21,:))+(cosd(120).*localinteraction(:,:,22,:))+(cosd(130).*localinteraction(:,:,23,:))+(cosd(140).*localinteraction(:,:,24,:))+(cosd(150).*localinteraction(:,:,25,:))+(cosd(160).*localinteraction(:,:,26,:))+(cosd(170).*localinteraction(:,:,27,:));
noverticaly12t(:,:,:) = (sind(0).*localinteraction(:,:,10,:))+(sind(10).*localinteraction(:,:,11,:))+(sind(20).*localinteraction(:,:,12,:))+(sind(30).*localinteraction(:,:,13,:))+(sind(40).*localinteraction(:,:,14,:))+(sind(50).*localinteraction(:,:,15,:))+(sind(60).*localinteraction(:,:,16,:))+(sind(70).*localinteraction(:,:,17,:))+(sind(80).*localinteraction(:,:,18,:))+(sind(90).*localinteraction(:,:,19,:))+(sind(-80).*localinteraction(:,:,20,:))+(sind(-70).*localinteraction(:,:,21,:))+(sind(-60).*localinteraction(:,:,22,:))+(sind(-50).*localinteraction(:,:,23,:))+(sind(-40).*localinteraction(:,:,24,:))+(sind(-30).*localinteraction(:,:,25,:))+(sind(-20).*localinteraction(:,:,26,:))+(sind(-10).*localinteraction(:,:,27,:));
nohorizontaly12t(:,:,:) = (sind(0).*localinteraction(:,:,10,:))+(sind(10).*localinteraction(:,:,11,:))+(sind(20).*localinteraction(:,:,12,:))+(sind(30).*localinteraction(:,:,13,:))+(sind(40).*localinteraction(:,:,14,:))+(sind(50).*localinteraction(:,:,15,:))+(sind(60).*localinteraction(:,:,16,:))+(sind(70).*localinteraction(:,:,17,:))+(sind(80).*localinteraction(:,:,18,:))+(sind(90).*localinteraction(:,:,19,:))+(sind(100).*localinteraction(:,:,20,:))+(sind(110).*localinteraction(:,:,21,:))+(sind(120).*localinteraction(:,:,22,:))+(sind(130).*localinteraction(:,:,23,:))+(sind(140).*localinteraction(:,:,24,:))+(sind(150).*localinteraction(:,:,25,:))+(sind(160).*localinteraction(:,:,26,:))+(sind(170).*localinteraction(:,:,27,:));

for u = 1:size(sortx)
noverticalxxt(sortx(u),sorty(u),:) = noverticalx12t(sortx(u),sorty(u),:);
noverticalyyt(sortx(u),sorty(u),:) = noverticaly12t(sortx(u),sorty(u),:);
nohorizontalxxt(sortx(u),sorty(u),:) = nohorizontalx12t(sortx(u),sorty(u),:);
nohorizontalyyt(sortx(u),sorty(u),:) = nohorizontaly12t(sortx(u),sorty(u),:);
end

for i = 1:1:neuroncol % this is the selection mechanism for orientation perception, for direction component, you need to also factor temporal component here
        for n = 1:1:neuronrow
            % get the beginning and ending column for this square
            endcol = i*squwidth;
            begcol = endcol - (squwidth-1);
            % get the beginning and ending row for this square
            endrow = n*squlength;
            begrow = endrow-(squlength-1);
            % get the data for that square
            readt = abs(sum(verticalinaccurate12t(begrow:endrow,begcol:endcol,:)));
            read1t = abs(sum(horizontalinaccurate12t(begrow:endrow,begcol:endcol,:)));
            
            if readt > read1t
                selectedx12t(n,i,:) = noverticalxxt(n,i,:);
                selectedy12t(n,i,:) = noverticalyyt(n,i,:);
            elseif readt < read1t
                selectedx12t(n,i,:) = nohorizontalxxt(n,i,:);
                selectedy12t(n,i,:) = nohorizontalyyt(n,i,:);
            end
vectorangles(n,i,:) = atan2d(selectedy12t(n,i,:),selectedx12t(n,i,:));
        end
end

roots = ceil(size(sortx,1)^0.5);
half = ceil(size(sortx,1)*0.5); % % for poggendorff illusion
root = ceil(half^0.5);
angleline = zeros(size(sortx,1),1);

for s= 2:size(sortx,1)-1
    angleline(s,:)=vectorangles(sortx(s),sorty(s),time);
    disp(angleline)
    
            % let's plot it
%            subplot(numlength,numwidth,(n-1)*numwidth+i)
%            plot(booksqueeze)
%            title([num2str(i) ' ' num2str(n)])
%             hypercolumnboxv1(s,:)=hypercolumnv1;
%              subplot(root,root,s); 
%              plot(subplotx,hypercolumnv1); title(['x' num2str(sortx(s)) ' '  'y' num2str(sorty(s))])
%              hold on; 
%              plot([1 1]*(90-localtheta), ylim, '--k')
%              hold off;
         output = struct;
    output.angle = angleline;
    
    save(['angles' 'interactionfieldsize' num2str(xsigmae) 'eiratio' num2str(eiratio) '.mat'],'-struct','output')        
             
end

% for s= 2:half
%     angleline(s,:)=vectorangles(sortx(s),sorty(s),time);
%     disp(angleline)
%     
%             % let's plot it
% %            subplot(numlength,numwidth,(n-1)*numwidth+i)
% %            plot(booksqueeze)
% %            title([num2str(i) ' ' num2str(n)])
% %             hypercolumnboxv1(s,:)=hypercolumnv1;
% %              subplot(root,root,s); 
% %              plot(subplotx,hypercolumnv1); title(['x' num2str(sortx(s)) ' '  'y' num2str(sorty(s))])
% %              hold on; 
% %              plot([1 1]*(90-localtheta), ylim, '--k')
% %              hold off;
%          output = struct;
%     output.angle = angleline;
%     
%     save(['angles' 'interactionfieldsize' num2str(xsigmae) 'eiratio' num2str(eiratio) '.mat'],'-struct','output')        
%              
% end
