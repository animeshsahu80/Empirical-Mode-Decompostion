
% EMD:  Emprical mode decomposition
% imf = emd(x)
% x   - input signal (must be a column vector)

function imf = emd(x);
fs = 1000;
t = 0:2/fs:2-1/fs;
x=0.8*sin(2*pi*50*t) + cos(2*pi*80*t)+ 0.6*sin(2*pi*25*t) + 0.4*cos(2*pi*10*t) +0.3*cos(2*pi*3*t);

c = x; % copy of the input signal 
N = length(x);

% loop to decompose the input signal into successive IMF

imf = []; % Matrix which will contain the successive IMF, and the residue

while (1) % the stop criterion is tested at the end of the loop
   
   
   % inner loop to find each imf
   
   h = c; % at the beginning of the sifting process, h is the signal
   SD = 1; % Standard deviation which will be used to stop the sifting process
   
   while SD > 0.3
      % while the standard deviation is higher than 0.3 
      
      % find local max/min points
      d = diff(h); % approximate derivative
      maxmin = []; % to store the optima 
      for i=1:N-2
         if d(i)==0                        % we are on a zero
            maxmin = [maxmin, i];
         elseif sign(d(i))~=sign(d(i+1))   % we are straddling a zero so
            maxmin = [maxmin, i+1];        % define zero as at i+1 
         end
      end
      
      if size(maxmin,2) < 2 % then it is the residue
         break
      end
      %.0000 + 1.0000i;
      % divide maxmin into maxes and mins
      if maxmin(1)>maxmin(2)              % first one is a max not a min
         maxes = maxmin(1:2:length(maxmin));
         mins  = maxmin(2:2:length(maxmin));
      else                                % is the other way around
         maxes = maxmin(2:2:length(maxmin));
         mins  = maxmin(1:2:length(maxmin));
      end
      
      % make endpoints both maxes and mins
      maxes = [1 maxes N];
      mins  = [1 mins  N];   
      
      % spline interpolate to get max and min envelopes; form imf
      maxenv = spline(maxes,h(maxes),1:N);
      minenv = spline(mins, h(mins),1:N);
      
      m = (maxenv + minenv)/2; % mean of max and min enveloppes
      prevh = h; % copy of the previous value of h before modifying it
      h = h - m; % substract mean to h
      
      % calculate standard deviation
      eps = 0.0000001; % to avoid zero values
      SD = sum ( ((prevh - h).^2) ./ (prevh.^2 + eps) );
      
   end
   
   imf = [imf; h]; % store the extracted IMF in the matrix imf
   % if size(maxmin,2)<2, then h is the residue
   
   % stop criterion of the algo.
   if size(maxmin,2) < 2
      break
   end
  
   c = c - h; % substract the extracted IMF from the signal
   
end
% subplot(4,length(imf)/4,1);
% plot(t,x);
% display(imf);
% for i=1:length(imf)
%     subplot(4,length(imf)/4,i+1);
%     plot(imf(i,1:1000));
% end
subplot(4,2,1);
plot(t,x);         %input signal
subplot(4,2,2);
plot(imf(1,1:1000));    %imf1
subplot(4,2,3);
plot(imf(2,1:1000));    %imf2
subplot(4,2,4);
plot(imf(3,1:1000));      %imf3
subplot(4,2,5);
plot(imf(4,1:1000));       %imf4
subplot(4,2,6);
plot(imf(5,1:1000));        %imf5(there are total 10 imf's and 1 residual function
% % xr=hilbert(imf(1,1:1000));
% % subplot(4,2,7);
% % plot(xr);

z = hilbert(imf(1,1:1000));     %Instantaneous frequency for 1st imf
instfreq = fs/(2*pi)*diff(unwrap(angle(z)));
%display(length(t));
display(length(instfreq));
subplot(4,2,7);
plot(t(2:1000),instfreq);
xlabel('Time')
ylabel('Hz')
grid on
title('Instantaneous Frequency')
z1 = hilbert(imf(2,1:1000));      %instantaneous frequency for 2nd imf
instfreq1 = fs/(2*pi)*diff(unwrap(angle(z1)));
subplot(4,2,8);
plot(t(2:1000),instfreq1);
xlabel('Time')
ylabel('Hz')
grid on
title('Instantaneous Frequency')
                                %by repeating above method we can find
                                %instantanous frequency for rest of the
                                %imsf's
return
