function [Theta,Phi,MapPhi,MapTheta]= beamforming(signalFreq, Pression, Position, c,npSamples)
% BEAMFORMING Executes beamforming method
DM = 0.0;
Phi = 0.0;
Theta = 0.0;
B = 0;
m = 0;
[x,n] = size(Pression);
Delta = zeros(1,4);
time = signalFreq^-1;
MapPhi = zeros(91, n+npSamples);
MapTheta = zeros(361, n+npSamples);
    for phi = 0:90,
        for theta = 0:360,
            for i = 1:4,
                dm = cosd(theta)*cosd(phi)*Position(i,1);
                dm = dm + cosd(phi)*sind(theta)*Position(i,2);
                dm = dm + sind(phi)*Position(i,3);
                Delta(i) = dm;
                if dm > DM
                    DM = dm;
                end
            end
            m = ceil(DM * npSamples );%/ signalFreq)
            Z = zeros(1,m);
            PressionAux = zeros(4,m+n);
            for i = 1:4,
                PressionAux(i,:) = cat(2,Z,Pression(i,:));
                if Delta(i) < 0
                    x = ceil(-1*Delta(i)*npSamples);%/signalFreq) %ver isso pois pode dar que sai do array
                    PressionAux(i,x+1:n+x) = PressionAux(i,m+1:m+n);
                    PressionAux(1,n+x+1:n+m) = zeros(1,m-x);
                else
                    x =  ceil(Delta(i)*npSamples);%/signalFreq); %var a possibilidade de x > 2
                    PressionAux(i,m-x+1:n+m-x) = PressionAux(i,m+1:n+m);
                    PressionAux(i,n+1+m-x:m+n) = zeros(1,x);
                end
                if i == 1 %tirar esse if
                    b = PressionAux(i,:);
                else
                    b = plus(b,PressionAux(i,:));
                end
            end  
            b(:) = b(:)/4; %normalizing the sum
            if max(b) > B
                B = max(b);
                Phi = phi;
                Theta = theta;
            end
            MapPhi(phi+1,1:m+n) = b;
            MapTheta(theta+1,1:m+n) = b;
        end   
    end
end