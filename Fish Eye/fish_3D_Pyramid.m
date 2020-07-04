% ������׶���ӽ�
% ˫���Բ�ֵ�Ż�
clc;
clear;
I0=imread('001.jpg');
I0=uint8(rgb2gray(I0));  
[heigh,width]=size(I0);
r=floor((min(heigh,width))/2);   %ԭͼ��뾶
oriwidth=2*r;                    %ԭͼ�����
oriheigh=2*r;
nwidth=4*r;                      %����У����ͼ���С
nheigh=4*r;
%ewidth1=round(0.707*r);
%ewidth2=round(1.293*r);
%eheigh1=round(0.707*r);
%eheigh2=round(1.293*r);
I1=uint8(zeros(nheigh,nwidth));
p1=r;
m1=r;
n1=r;
xc=r;
yc=r;
u0=r;
v0=r;

% �����ȡ��������k1(����)��k2�����£�
k1=round(1.65*r);
k2=round(0.4*r);

%-------------����У��-------%

for i=1:2*r  
    for j=k1:k1+2*r
%         m=2*r-j;
%         n=r-i;
        m=j-2*r;
        n=i-r;
        u=(r*n)/(sqrt(0.5*m^2+n^2+(2*r^2-2*r*m+0.5*m^2)));
        v=(0.7*r*m)/(sqrt(0.5*m^2+n^2+(2*r^2-2*r*m+0.5*m^2)));
        xft=round(u+r);                        
        yft=round(v+r);
          if (xft<=0)||(xft>oriwidth)||(yft<=0)||(yft>oriheigh)
            I1(i,j,:)=0;
          % ˫���Բ�ֵ��ϸ���Ż�
          elseif (xft >0) && (xft <oriwidth) && (yft>0) && (yft <oriheigh)            
              xx= u+u0;
              yy=v+v0;
              a = double(uint16(xx)); 
              b = double(uint16(yy)); 
              x11 = double(I0(a,b));                     % x11 <- origImg(a,b)
              x12 = double(I0(a,b+1));                   % x12 <- origImg(a,b+1)
              x21 = double(I0(a+1,b));                   % x21 <- origImg(a+1,b)
              x22 = double(I0(a+1,b+1));                 % x22 <- origImg(a+1,b+1)
              I1(i,j+2*r-k1) = uint8( (b+1-yy) * ((xx-a)*x21 + (a+1-xx)*x11) + (yy-b) * ((xx-a)*x22 +(a+1-xx) * x12) ); %˫���Բ�ֵ
          else                                                      
              I1(i,j+2*r-k1) = I0(xft,yft);     % else,apply nearest interpolate
          end      
    end  
end

%-------------����У��-------%
for i=1:2*r  
    for j=2*r-k1:4*r-k1
        m=2*r-j;
        n=i-r;
        u=(r*n)/(sqrt(0.5*m^2+n^2+(2*r^2-2*r*m+0.5*m^2)));
        v=-(0.7*r*m)/(sqrt(0.5*m^2+n^2+(2*r^2-2*r*m+0.5*m^2)));
        xft=round(u+r);                        
        yft=round(v+r);
        if (xft<=0)||(xft>oriwidth)||(yft<=0)||(yft>oriheigh)
            I1(i,j,:)=0;
        elseif (xft >0) && (xft <oriwidth) && (yft>0) && (yft <oriheigh)
            xx= u+u0;
            yy=v+v0;
            a = double(uint16(xx)); 
            b = double(uint16(yy)); 
            x11 = double(I0(a,b));                     % x11 <- origImg(a,b)
            x12 = double(I0(a,b+1));                   % x12 <- origImg(a,b+1)
            x21 = double(I0(a+1,b));                   % x21 <- origImg(a+1,b)
            x22 = double(I0(a+1,b+1));                 % x22 <- origImg(a+1,b+1)
            I1(i,j+k1-2*r+1) = uint8( (b+1-yy) * ((xx-a)*x21 + (a+1-xx)*x11) + (yy-b) * ((xx-a)*x22 +(a+1-xx) * x12) ); %˫���Բ�ֵ
        else
            I1(i,j+k1-2*r+1) = I0(xft,yft);     % else,apply nearest interpolate        
        end     
    end  
end

%-------------����У��-------%
for i=k2:k2+2*r
    for j=1:2*r
        m=j-r;
        n=2*r-i;
        u=-(0.7*r*n)/(sqrt(0.5*n^2+m^2+(2*r^2-2*r*n+0.5*n^2)));
        v=(r*m)/(sqrt(0.5*n^2+m^2+(2*r^2-2*r*n+0.5*n^2)));
        xft=round(u+r);                        
        yft=round(v+r);
        if (xft<=0)||(xft>oriwidth)||(yft<=0)||(yft>oriheigh)
            I1(i,j,:)=0;
        elseif (xft >0) && (xft <oriwidth) && (yft>0) && (yft <oriheigh)
            xx= u+u0;
            yy=v+v0;
            a = double(uint16(xx)); 
            b = double(uint16(yy)); 
            x11 = double(I0(a,b));                     % x11 <- origImg(a,b)
            x12 = double(I0(a,b+1));                   % x12 <- origImg(a,b+1)
            x21 = double(I0(a+1,b));                   % x21 <- origImg(a+1,b)
            x22 = double(I0(a+1,b+1));                 % x22 <- origImg(a+1,b+1)
            I1(i+2*r-k2+1,j) = uint8( (b+1-yy) * ((xx-a)*x21 + (a+1-xx)*x11) + (yy-b) * ((xx-a)*x22 +(a+1-xx) * x12) ); %˫���Բ�ֵ
        else    
            I1(i+2*r-k2+1,j) = I0(xft,yft);     % else,apply nearest interpolate    
        end      
    end  
end

%-------------����У��-------%
for i=2*r-k2:4*r-k2
    for j=1:2*r
        m=r-j;
        n=i-2*r;
        u=(0.7*r*n)/(sqrt(0.5*n^2+m^2+(2*r^2-2*r*n+0.5*n^2)));
        v=-(r*m)/(sqrt(0.5*n^2+m^2+(2*r^2-2*r*n+0.5*n^2)));
        xft=round(u+r);                        
        yft=round(v+r);
        if (xft<=0)||(xft>oriwidth)||(yft<=0)||(yft>oriheigh)
            I1(i,j,:)=0;
        elseif (xft >0) && (xft <oriwidth) && (yft>0) && (yft <oriheigh)
            xx= u+u0;
            yy=v+v0;
            a = double(uint16(xx)); 
            b = double(uint16(yy)); 
            x11 = double(I0(a,b));                     % x11 <- origImg(a,b)
            x12 = double(I0(a,b+1));                   % x12 <- origImg(a,b+1)
            x21 = double(I0(a+1,b));                   % x21 <- origImg(a+1,b)
            x22 = double(I0(a+1,b+1));                 % x22 <- origImg(a+1,b+1)
            I1(i+k2+1,j+2*r) = uint8( (b+1-yy) * ((xx-a)*x21 + (a+1-xx)*x11) + (yy-b) * ((xx-a)*x22 +(a+1-xx) * x12) ); %˫���Բ�ֵ
        else               
                I1(i+k2+1,j+2*r) = I0(xft,yft);     % else,apply nearest interpolate        
        end      
    end  
end

%figure,imshow(I1);
figure;
subplot(1,2,1),imshow(I0);
subplot(1,2,2),imshow(I1);