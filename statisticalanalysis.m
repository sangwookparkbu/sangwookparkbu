 
fileName = ['angles' 'interactionfieldsize' num2str(1) 'eiratio' num2str(1.1) '.mat'];
   
    loading1 = load(fileName);
    loading1.angle(1:5,:) = [];
    loading1.angle(end-4:end,:) = [];
   
fileName2 = ['angles' 'interactionfieldsize' num2str(1.2) 'eiratio' num2str(1.1) '.mat'];
    
    loading2 = load(fileName2);
    loading2.angle(1:5,:) = [];
    loading2.angle(end-4:end,:) = [];
    
fileName3 = ['angles' 'interactionfieldsize' num2str(1.4) 'eiratio' num2str(1.1) '.mat'];
    
    loading3 = load(fileName3);
    loading3.angle(1:5,:) = [];
    loading3.angle(end-4:end,:) = [];
 
    
fileName4 = ['angles' 'interactionfieldsize' num2str(1.6) 'eiratio' num2str(1.1) '.mat'];
   
    loading4 = load(fileName4);
    loading4.angle(1:5,:) = [];
    loading4.angle(end-4:end,:) = [];
    
fileName5 = ['angles' 'interactionfieldsize' num2str(1.8) 'eiratio' num2str(1.1) '.mat'];
   
    loading5 = load(fileName5);
    loading5.angle(1:5,:) = [];
    loading5.angle(end-4:end,:) = [];

fileName6 = ['angles' 'interactionfieldsize' num2str(2) 'eiratio' num2str(1.1) '.mat'];
   
    loading6 = load(fileName6);
    loading6.angle(1:5,:) = [];
    loading6.angle(end-4:end,:) = [];

fileName7 = ['angles' 'interactionfieldsize' num2str(2.2) 'eiratio' num2str(1.1) '.mat'];
   
    loading7 = load(fileName7);
    loading7.angle(1:5,:) = [];
    loading7.angle(end-4:end,:) = [];

fileName8 = ['angles' 'interactionfieldsize' num2str(2.4) 'eiratio' num2str(1.1) '.mat'];
   
    loading8 = load(fileName8);
    loading8.angle(1:5,:) = [];
    loading8.angle(end-4:end,:) = [];
    
fileName9 = ['angles' 'interactionfieldsize' num2str(2.6) 'eiratio' num2str(1.1) '.mat'];
   
    loading9 = load(fileName9);
    loading9.angle(1:5,:) = [];
    loading9.angle(end-4:end,:) = [];

fileName10 = ['angles' 'interactionfieldsize' num2str(2.8) 'eiratio' num2str(1.1) '.mat'];
   
    loading10 = load(fileName10);
    loading10.angle(1:5,:) = [];
    loading10.angle(end-4:end,:) = [];

fileName11 = ['angles' 'interactionfieldsize' num2str(3) 'eiratio' num2str(1.1) '.mat'];
   
    loading11 = load(fileName11);
    loading11.angle(1:5,:) = [];
    loading11.angle(end-4:end,:) = [];
    
fileName12 = ['angles' 'interactionfieldsize' num2str(3.2) 'eiratio' num2str(1.1) '.mat'];
   
    loading12 = load(fileName12);
    loading12.angle(1:5,:) = [];
    loading12.angle(end-4:end,:) = [];
    
    
    lines = zeros(size(loading1.angle,1),12);
    
    lines(:,1) = loading1.angle;
    lines(:,2) = loading2.angle;
    lines(:,3) = loading3.angle;
    lines(:,4)= loading4.angle;
    lines(:,5)=loading5.angle;
    lines(:,6) = loading6.angle;
    lines(:,7) = loading7.angle;
    lines(:,8) = loading8.angle;
    lines(:,9) = loading9.angle;
    lines(:,10) = loading10.angle;
    lines(:,11) = loading11.angle;
    lines(:,12) = loading12.angle;
    % lines(:,4) = loading6.angle;
%     lines = lines - 70;

     plot(lines(:,1)); hold on; plot(lines(:,2)); hold on; plot(lines(:,3)); hold on; plot(lines(:,4)); hold on; plot(lines(:,5)); hold on; plot(lines(:,6)); hold on; plot(lines(:,7)); hold on; plot(lines(:,8)); hold on; plot(lines(:,9)); hold on; plot(lines(:,10)); hold on; plot(lines(:,11)); hold on; plot(lines(:,12)); legend('Interaction field size 1','Interaction field size 2','Interaction field size 3'); xlabel('Spatial Location'); ylabel('Angle Response');
     
     lines2 = lines - 70; % subtracting the orientation of the tilted line
     
     plot(lines2(:,1)); hold on; plot(lines2(:,2)); hold on; plot(lines2(:,3)); hold on; plot(lines2(:,4)); hold on; plot(lines2(:,5)); hold on; plot(lines2(:,6)); hold on; plot(lines2(:,7)); hold on; plot(lines2(:,8)); hold on; plot(lines2(:,9)); hold on; plot(lines2(:,10)); hold on; plot(lines2(:,11)); hold on; plot(lines2(:,12)); legend('Interaction field size 1','Interaction field size 2','Interaction field size 3'); xlabel('Spatial Location'); ylabel('Angle Response');
     
%      threshold = -5; % threshold for magnitude of illusion
%      
%      count1 = sum(lines2(:,1) < threshold);
%      count2 = sum(lines2(:,2) < threshold);
%      count3 = sum(lines2(:,3) < threshold);
%      
%      count11 = find(lines2(:,1) < threshold);
%      count11(:,:) = 1;
%      count22 = find(lines2(:,2) < threshold);
%      count22(:,:) = 2;
%      count33 = find(lines2(:,3) < threshold);
%      count33(:,:) = 3;
     
     firstrepulsion = zeros(size(lines2,1),1);
     firstattraction = zeros(size(lines2,1),1);
     secondrepulsion = zeros(size(lines2,1),1);
     secondattraction = zeros(size(lines2,1),1);
     thirdrepulsion = zeros(size(lines2,1),1);
     thirdattraction = zeros(size(lines2,1),1);
      fourthrepulsion = zeros(size(lines2,1),1);
      fourthattraction = zeros(size(lines2,1),1);
      fifthrepulsion = zeros(size(lines2,1),1);
      fifthattraction = zeros(size(lines2,1),1);
      sixthrepulsion = zeros(size(lines2,1),1);
      sixthattraction = zeros(size(lines2,1),1);
      seventhrepulsion = zeros(size(lines2,1),1);
      seventhattraction = zeros(size(lines2,1),1);
      eighthrepulsion = zeros(size(lines2,1),1);
      eighthattraction = zeros(size(lines2,1),1);
      ninethrepulsion = zeros(size(lines2,1),1);
      ninethattraction = zeros(size(lines2,1),1);
      tenthrepulsion = zeros(size(lines2,1),1);
      tenthattraction = zeros(size(lines2,1),1);
      eleventhrepulsion = zeros(size(lines2,1),1);
      eleventhattraction =  zeros(size(lines2,1),1);
      twelvethrepulsion = zeros(size(lines2,1),1);
      twelvethattraction = zeros(size(lines2,1),1);
     
     for a = 1:size(lines2,1)
         if lines2(a,1)<0
             firstrepulsion(a,1) = lines2(a,1);
         elseif lines2(a,1)>0
             firstattraction(a,1) = lines2(a,1);
         end
         
         if lines2(a,2)<0
             secondrepulsion(a,1) = lines2(a,2);
         elseif lines2(a,2)>0
             secondattraction(a,1) = lines2(a,2);
         end
         
         if lines2(a,3)<0
             thirdrepulsion(a,1) = lines2(a,3);
         elseif lines2(a,3)>0
             thirdattraction(a,1) = lines2(a,3);
         end
         
         if lines2(a,4)<0
             fourthrepulsion(a,1) = lines2(a,4);
         elseif lines2(a,4)>0
             fourthattraction(a,1) = lines2(a,4);
         end
         
         if lines2(a,5)<0
             fifthrepulsion(a,1) = lines2(a,5);
         elseif lines2(a,5)>0
             fifthattraction(a,1) = lines2(a,5);
         end
         
         if lines2(a,6)<0
             sixthrepulsion(a,1) = lines2(a,6);
         elseif lines2(a,6)>0
             sixthattraction(a,1) = lines2(a,6);
         end
         
         if lines2(a,7)<0
             seventhrepulsion(a,1) = lines2(a,7);
         elseif lines2(a,7)>0
             seventhattraction(a,1) = lines2(a,7);
         end
         
         if lines2(a,8)<0
             eighthrepulsion(a,1) = lines2(a,8);
         elseif lines2(a,8)>0
             eighthattraction(a,1) = lines2(a,8);
         end
         
         if lines2(a,9)<0
             ninethrepulsion(a,1) = lines2(a,9);
         elseif lines2(a,9)>0
             ninethattraction(a,1) = lines2(a,9);
         end
         
         if lines2(a,10)<0
             tenthrepulsion(a,1) = lines2(a,10);
         elseif lines2(a,10)>0
             tenthattraction(a,1) = lines2(a,10);
         end
         
         if lines2(a,11)<0
             eleventhrepulsion(a,1) = lines2(a,11);
         elseif lines2(a,11)>0
             eleventhattraction(a,1) = lines2(a,11);
         end
         
          if lines2(a,12)<0
             twelvethrepulsion(a,1) = lines2(a,12);
         elseif lines2(a,12)>0
             twelvethattraction(a,1) = lines2(a,12);
         end
     end
             
     firstrepulsionmagnitude = abs(sum(firstrepulsion));
     firstattractionmagnitude = abs(sum(firstattraction));
     secondrepulsionmagnitude = abs(sum(secondrepulsion));
     secondattractionmagnitude = abs(sum(secondattraction));
     thirdrepulsionmagnitude = abs(sum(thirdrepulsion));
     thirdattractionmagnitude = abs(sum(thirdattraction));
      fourthrepulsionmagnitude = abs(sum(fourthrepulsion));
      fourthattractionmagnitude = abs(sum(fourthattraction));
      fifthrepulsionmagnitude = abs(sum(fifthrepulsion));
      fifthattractionmagnitude = abs(sum(fifthattraction));
      sixthrepulsionmagnitude = abs(sum(sixthrepulsion));
      sixthattractionmagnitude = abs(sum(sixthattraction));
      seventhrepulsionmagnitude = abs(sum(seventhrepulsion));
      seventhattractionmagnitude = abs(sum(seventhattraction));
      eighthrepulsionmagnitude = abs(sum(eighthrepulsion));
      eighthattractionmagnitude = abs(sum(eighthattraction));
      ninethrepulsionmagnitude = abs(sum(ninethrepulsion));
      ninethattractionmagnitude = abs(sum(ninethattraction));
      tenthrepulsionmagnitude = abs(sum(tenthrepulsion));
      tenthattractionmagnitude = abs(sum(tenthattraction));
      eleventhrepulsionmagnitude = abs(sum(eleventhrepulsion));
      eleventhattractionmagnitude = abs(sum(eleventhattraction));
      twelvethrepulsionmagnitude = abs(sum(twelvethrepulsion));
      twelvethattractionmagnitude = abs(sum(twelvethattraction));
      
     ratio1 = firstrepulsionmagnitude/firstattractionmagnitude;
     ratio2 = secondrepulsionmagnitude/secondattractionmagnitude;
     ratio3 = thirdrepulsionmagnitude/thirdattractionmagnitude;
     ratio4 = fourthrepulsionmagnitude/fourthattractionmagnitude;
     ratio5 = fifthrepulsionmagnitude/fifthattractionmagnitude;
     ratio6 = sixthrepulsionmagnitude/sixthattractionmagnitude;
     ratio7 = seventhrepulsionmagnitude/seventhattractionmagnitude;
     ratio8 = eighthrepulsionmagnitude/eighthattractionmagnitude;
     ratio9 = ninethrepulsionmagnitude/ninethattractionmagnitude;
     ratio10 = tenthrepulsionmagnitude/tenthattractionmagnitude;
     ratio11 = eleventhrepulsionmagnitude/eleventhattractionmagnitude;
     ratio12 = twelvethrepulsionmagnitude/twelvethattractionmagnitude;
     
%      ratio4 = fourthrepulsionmagnitude/fourthattractionmagnitude;
     
%      repulsionhisto1 = sum(firstrepulsion < 0);
%      repulsionhisto2 = sum(secondrepulsion < 0);
%      repulsionhisto3 = sum(thirdrepulsion < 0);
%      attractionhisto1 = sum(firstattraction > 0);
%      attractionhisto2 = sum(secondattraction > 0);
%      attractionhisto3 = sum(thirdattraction > 0);
     
     repulsionhisto1 = zeros(round(firstrepulsionmagnitude),1);
     repulsionhisto1(:,:) = 1;
     attractionhisto1 = zeros(round(firstattractionmagnitude),1);
     attractionhisto1(:,:) = 1;
     repulsionhisto2 = zeros(round(secondrepulsionmagnitude),1);
     repulsionhisto2(:,:) = 2;
     attractionhisto2 = zeros(round(secondattractionmagnitude),1);
     attractionhisto2(:,:) = 2;
     repulsionhisto3 = zeros(round(thirdrepulsionmagnitude),1);
     repulsionhisto3(:,:) = 3;
     attractionhisto3 = zeros(round(thirdattractionmagnitude),1);
     attractionhisto3(:,:) = 3;
     
     
%      histoplot = zeros(3,1);
% 
%      histoplot(1,1) = count1;
%      histoplot(2,1) = count2;
%      histoplot(3,1) = count3;
     repulsioncorrelationx = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2];  % size of interaction field from 1 to 3
     repulsioncorrelationy = [firstrepulsionmagnitude secondrepulsionmagnitude thirdrepulsionmagnitude fourthrepulsionmagnitude fifthrepulsionmagnitude sixthrepulsionmagnitude seventhrepulsionmagnitude eighthrepulsionmagnitude ninethrepulsionmagnitude tenthrepulsionmagnitude eleventhrepulsionmagnitude twelvethrepulsionmagnitude]; % magnitude of repulsion 
     repulsionslopecal = polyfit(repulsioncorrelationx, repulsioncorrelationy, 1);
     repulsionslope = repulsionslopecal(1);
     repulsiongraph=fitlm(repulsioncorrelationx,repulsioncorrelationy);             % a linear model
     [repulsioncorrelationcoeffecient, repulsionpvalue] = corrcoef(repulsioncorrelationx,repulsioncorrelationy);
     
     attractioncorrelationx = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2];  % size of interaction field from 1 to 3
     attractioncorrelationy = [firstattractionmagnitude secondattractionmagnitude thirdattractionmagnitude fourthattractionmagnitude fifthattractionmagnitude sixthattractionmagnitude seventhattractionmagnitude eighthattractionmagnitude ninethattractionmagnitude tenthattractionmagnitude eleventhattractionmagnitude twelvethattractionmagnitude]; % magnitude of repulsion 
     attractionslopecal = polyfit(attractioncorrelationx, attractioncorrelationy, 1);
     attractionslope = attractionslopecal(1);
     attractiongraph=fitlm(attractioncorrelationx,attractioncorrelationy);             % a linear model
     attractioncorrelationcoeffecient = corr2(attractioncorrelationx,attractioncorrelationy);
%      histoplot1 = zeros(3,1);
%      histoplot1(1,1) = count11;
%      histoplot1(2,1) = count22;
%      histoplot1(3,1) = count33;
ratioanalysis = [ratio1 ratio2 ratio3 ratio4 ratio5 ratio6 ratio7 ratio8 ratio9 ratio10 ratio11 ratio12];
ratiograph=fitlm(repulsioncorrelationx,ratioanalysis); 
ratiobar = bar(repulsioncorrelationx, ratioanalysis); hold on; plot(ratiograph); hold on; xtips1 = ratiobar(1).XEndPoints;
ytips1 = ratiobar(1).YEndPoints;
labels1 = string(ratiobar(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom');

            % a linear model
     [ratiocorrelationcoeffecient, ratiopvalue] = corrcoef(repulsioncorrelationx,ratioanalysis);
     
     
 

bary = [repulsioncorrelationy; attractioncorrelationy];
Illusionmagnitudebar = bar(attractioncorrelationx,bary); hold on; xtips1 = Illusionmagnitudebar(1).XEndPoints;
ytips1 = Illusionmagnitudebar(1).YEndPoints;
labels1 = string(Illusionmagnitudebar(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom'); xtips2 = Illusionmagnitudebar(2).XEndPoints;
ytips2 = Illusionmagnitudebar(2).YEndPoints;
labels2 = string(Illusionmagnitudebar(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
      'VerticalAlignment','bottom'); hold on; p = plot(repulsiongraph); ylim([0 thirdrepulsionmagnitude*1.5]); xlabel('Interaction Field Size'); ylabel('Magnitude of Distortion'); title('Size of Interaction Field vs Magnitude of Distortion');
   % p(end-1,1).Visible='off'; p(end,1).Visible='off'  
    % histogram(repulsionhisto1,'BinEdges',[0.75 1]); hold on; histogram(attractionhisto1,'BinEdges',[1 1.25]); hold on; histogram(repulsionhisto2,'BinEdges',[1.75 2]); hold on; histogram(attractionhisto2,'BinEdges',[2 2.25]); hold on; histogram(repulsionhisto3,'BinEdges',[2.75 3]); hold on; histogram(attractionhisto3,'BinEdges',[3 3.25]); hold on; p = plot(repulsiongraph); p(end-1,1).Visible='off'; p(end,1).Visible='off'; ylim([0 thirdrepulsionmagnitude*1.5]); xlabel('Interaction Field Size'); ylabel('Magnitude of Distortion'); title('Size of Interaction Field vs Magnitude of Distortion');
        

%      
%      lines3(:,1) = loading2.angle;
%      lines3(:,2) = loading4.angle;
%      lines3(:,3) = loading5.angle;
%      
%      plot(lines3(:,1)); hold on; plot(lines3(:,2)); hold on; plot(lines3(:,3)); legend('eiratio1/1.1','eiratio1/1.6','eiratio1/2.1'); xlabel('spatial location'); ylabel('angle response to tilted line');
%      
%      [pval,tbll,statss] = kruskalwallis(lines3);
%     cc = multcompare(statss);
    
%     [p,tbl,stats] = kruskalwallis(lines);
%     c = multcompare(stats);
%     
%     predictor = [1,2,3];
%     output = [3, 3 ,3];
%     
%     y = zeros(size(lines,1),3);
%     mvregress
%     plot(predictor,lines)
%     
%     fine = zeros(size(loading1.angle,1),2);
%     fine(:,1) = lines(:,1);
%     fine(:,2) = 70;
%     
   
    