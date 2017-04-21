function [ amplitude,phase ] = nanospectrum(frequency,e1,e2)
%this function calculates the amplitude and phase of the far field 
%scattering of the apertureless near field microscopy. 
%   frequency is the column vector of sampled points in frequency range. 
%   e1 is the column vector of the sample dielectric permittivity e1.
%   e2 is the column vector of the sample dielectric permittivity e2.
%   amplitude is the column vector of the scattered wave of the tip, the
%   demodulation order is determined by the parameter within the codes.
%   phase is the column vector of the scattered wave of the tip, the
%   demodulation order is determined by the parameter within the codes.
n=size(frequency,1);
n1=500;                           %n1 is the sampling rate in time period
s=zeros(n,n1);
et=-1300+1i*760;                  %et is the dielectric function of the tip
a=15*10^(-9);                     %a is the radius of the tip
A0=65*10^(-9);                    %A0 is amplitude of the tip oscillation
alpha=4*pi*a^3*(et-1)/(et+2);
beta=(e1+1i.*e2-1)./(e1+1i.*e2+1);
f=2*pi*237579;                    %f is the oscillation frequency of the 
                                  %tip.
m=2;                              %m is the demodulation order of the 
                                  %signals.
inc=pi/4;                             %inc is the incident angle of the 
                                      %focused beam on the tip, the angle
                                      %between the tip and incident beam.
                                      %inc is in radians.



nr=sqrt(0.5*sqrt(e1.^2+e2.^2)+0.5*e1); %nr is Re(refractive index)
ki=sqrt(0.5*sqrt(e1.^2+e2.^2)-0.5*e1); %ki is Im(refractive index)

rf=((nr+1i.*ki).^2*cos(inc)-sqrt((nr+1i.*ki).^2-sin(inc)^2))./ ...
   ((nr+1i.*ki).^2*cos(inc)+sqrt((nr+1i.*ki).^2-sin(inc)^2));                             
%this equation calculates rf from el and e2.

%rf=rp.*exp(-1i.*phip*pi/180);

%rf is the reflection coefficient in far field, rp is the amplitude of the 
%fresnel reflectance parallel to the POI and phip is its phase. This takes
%the phase from WVASE, makes it positive, and converts it to radians so it
%will work in matlab.

p=1:n1;
t=2*pi/f/(n1-1)*(p-1);          
z=A0*(1+cos(f*t));              

for l=1:n1
    s(:,l)=alpha/(4*pi*a^3)./(1-alpha.*beta/...
        (16*pi*(a+z(l))^3))*exp(1i*m*f*t(l));
end


%Sm=trapz(t,s,2);
Sm=trapz(t,s,2).*(1+rf).^2;

amplitude=abs(Sm);
phase=angle(Sm);

subplot(2,2,1)
plot(frequency,amplitude)
subplot(2,2,2)
plot(frequency,phase)



%subplot(2,2,3)
%plot(frequency,rp)
%subplot(2,2,4)
%plot(frequency,-phip)

%once again the phase is made positive to adhere to the physics convention.
%(it is negative in the nebraska convention).

end