% scrip para despues de experimentoNOISEHF3
%clear all 
for ii=1:length(R3)
datos=R3{ii};
l=size(datos); 
l2=l(2);



for i=1:l2
    d1= cell2mat(datos(:,i));
     b1=[]; b2=[]; b3=[]; b4=[]; b5=[]; b6=[]; b7=[]; b8=[]; b9=[]; b10=[]; b11=[]; b12=[];
     %b13=[]; b14=[]; b15=[]; b16=[]; b17=[]; b18=[]; b19=[]; b20=[]; b21=[]; b22=[]; b23=[]; b24=[]; 
     %b25=[]; b26=[]; b27=[]; b28=[]; b29=[]; b30=[]; b31=[]; b32=[]; b33=[]; b34=[]; b35=[]; b36=[];
    
    for ii=1:length(d1)
        d2= d1(ii);

        if d2<=92.66
            b1(ii)=d2;
        elseif d2<=111.33
            b2(ii)=d2;
        elseif d2<= 130
            b3(ii)=d2;
        elseif d2<=148.66
            b4(ii)=d2;
        elseif d2<=167.33
            b5(ii)=d2;
        elseif d2<=186
            b6(ii)=d2;
        elseif d2<=204.666
            b7(ii)=d2;
        elseif d2<=223.33
            b8(ii)=d2;
        elseif d2<=242
            b9(ii)=d2;
        elseif d2<=260.666
            b10(ii)=d2;
        elseif d2<=279.33
            b11(ii)=d2;
        elseif d2<=298
            b12(ii)=d2;
        end 

    end

% B1(i)= sum(1<=b1); B2(i)= sum(1<=b2); B3(i)= sum(1<=b3); B4(i)= sum(1<=b4);
% B5(i)= sum(1<=b5); B6(i)= sum(1<=b6); B7(i)= sum(1<=b7); B8(i)= sum(1<=b8);
% B9(i)= sum(1<=b9); B10(i)= sum(1<=b10); B11(i)= sum(1<=b11); B12(i)= sum(1<=b12);

B1= sum(1<=b1); B2= sum(1<=b2); B3= sum(1<=b3); B4= sum(1<=b4);
B5= sum(1<=b5); B6= sum(1<=b6); B7= sum(1<=b7); B8= sum(1<=b8);
B9= sum(1<=b9); B10= sum(1<=b10); B11= sum(1<=b11); B12= sum(1<=b12);
R= [B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12];

R2(i,ii)= max(R);
end
end

figure
 plot(R2,'.')






'terminado :D'