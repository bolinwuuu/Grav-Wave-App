% Script to generate simple chirp waveform from within Matlab app designer

 % K. Riles - Summer 2022

% Many of the equations used in the simulation come from the following article:
% K. Riles, Prog. Part. Nucl. Phys. 68 (2013) 1-54 • e-Print: https://arxiv.org/abs/1209.0667
% (referred to below as "PPNP")

function status = run_chirp(Axes,Mass1, Mass2, RunMode)
  
  status = 0;
  fprintf('Entered run_chirp with\n Mass1 = %f\n Mass2 = %f\n RunMode = %d\n',Mass1,Mass2,RunMode);

  M1 = Mass1;
  M2 = Mass2;

  % Implied chirp mass (governs frequency and amplitude evolution)
  % (PPNP text right after Eqn 74)
  
  Mchirp = (M1*M2)^(3/5)/(M1+M2)^(1/5);

  % Physical constants
 
  G = 6.67e-11;
  c = 2.998e8;
  pc = 3.086e16;
  Msun = 2.0e30;

  % Compute Schwarzchild radii of stars

  R1 = (2*G*M1*Msun/c^2);
  R2 = (2*G*M2*Msun/c^2);

  % Frequency coefficient
  % (Based on PPNP Eqn 73)

  fcoeff = (1/(8*pi)) * (5^3)^(1/8) * (c^3/(G*Mchirp*Msun))^(5/8);

  % Amplitude coefficient (assume source at 15 Mpc)
  % (Based on PPNP Eqn 74)

  rMpc = 15.;
  r = rMpc * 1.e6 * pc;
  hcoeff = (1/r) * (5*G^5*(Mchirp*Msun)^5/c^11)^(1/4);

  % Amplitude rescaling parameter

  hscale = 1.0e21;

  % frequency (Hz) when signal enters detector band

  fbandlo = 30.; 

  % Compute time remaining to coalescence from entering band
  % (Based on PPNP Eqn 73)

  tau = (fcoeff/fbandlo)^(8/3);

  % Debugging summary

  fprintf('Starting chirp simulation with M1, M2, Mchirp = %.2f %.2f %.2f (Msun)\n',M1,M2,Mchirp);
  fprintf('--> Schwarzchild radii = %f %f m\n',R1,R2);
  fprintf('Distance to source r = %.1f Mpc\n',rMpc);
  fprintf('Detection band low frequency = %f Hz\n--> Time to coalescence = %f s\n',fbandlo,tau);

  % Sampling rate (Hz) - fixed at 48 kHz for mp4 output

  fsamp = 48000;   
  dt = 1/fsamp;

  % Length of time to simulate (round up to nearest thirtieth of an integer second and add a thirtieth)

  T = ceil(30*tau)/ 30. + min(1/10.,0.1*tau);

  % Create time sample container

  N = floor(fsamp*T); 
  t = [0:N-1].'*dt; 

  % Determine frequency (and then time) when Schwarzchild radii touch
  % (Use Kepler's 3rd law)
  % (Double orbital frequency to get GW frequency)

  ftouch = 2 * (1/(2*pi)) * (G*(M1+M2)*Msun/(R1+R2)^3)^(1/2);
  tautouch = (fcoeff/ftouch)^(8/3);

  fprintf('GW frequency when Schwarzchild radii touch: %.2f Hz\n--> Occurs %.6f seconds before point-mass coalescence\n',ftouch,tautouch);

  % Create frequency value vs time (up to last time sample before point-mass coalescence)
  % (Based on PPNP Eqn 73)

  freq = t*0;
  freq(t<tau) = fcoeff * (tau-t(t<tau)).^(-3/8);

  % Determine last valid frequency

  [maxfreq lastsample] = max(freq(freq<ftouch));
  freqlast = freq(lastsample);
  fprintf('Maximum frequency at time sample %d = %f Hz\n',lastsample,freqlast);

  % Fill frequency array with constant values for later times
  % for use with windowed audio waveform

  freq(lastsample+1:N) = freqlast;

  % Create amplitude value vs time (up to last time sample before touch)
  % (Based on PPNP Eqn 74)

  amp = t*0;
  amp(1:lastsample) = hcoeff * hscale * (tau-t(1:lastsample)).^(-1/4);

  % Generate strain signal in time domain

  h = t*0;

  % The following model is too simplistic:

  % h = amp .* sin(2*pi*freq.*t);

  % Instead, integrate phase numerically from frequency

  phi = t*0;
  phi = 2*pi*cumsum(freq)*dt;
  h = amp .* sin(phi);

  % Plot waveform value vs time

  if (RunMode==2)
    
    cla(Axes);
    plot(Axes,t, h);
    set(Axes,'YScale','linear')

    xlabel(Axes,'Time (s)');
    ylabel(Axes,'h(t) [arb. units]');
    %titlestr = sprintf('Inspiral waveform in time domain (%.1f,%.1f)',M1,M2);
    %title(titlestr);
    xlim(Axes,[0 T]);
    ylim(Axes,[-max(abs(h))*1.1 max(abs(h))*1.1]);
    %fname_png = sprintf('Chirp_waveform_%d_%d.png',floor(M1),floor(M2));
    %print('-dpng',fname_png);
    return;
    
  end
  
  % Plot waveform frequency vs time
  
  if (RunMode==3)
    
    cla(Axes);
    semilogy(Axes,t, freq);
    xlabel(Axes,'Time (s)');
    ylabel(Axes,'Frequency (Hz)');
    %titlestr = sprintf('Inspiral frequency vs time (%.1f,%.1f)',M1,M2);
    %title(titlestr);
    xlim(Axes,[0 T]);
    ylim(Axes,[fbandlo/3 ftouch*2])
    %fname_png = sprintf('Chirp_frequency_%d_%d.png',floor(M1),floor(M2)); 
    %print('-dpng',fname_png);
    return;

  end
  
% Plot spectrogram
  
  if (RunMode==4)
    
    cla(Axes);
    Nfft = round(floor(N/2)/200)*200;
    fprintf('Creating spectrogram with Nfft = %d\n',Nfft);
    [S,Fspect,Tspect] = spectrogram(h,Nfft/20,0,Nfft/10,fsamp);
    %spectrogram(h,Nfft/100,0,Nfft/100,fsamp,'yaxis');
    surf(Axes,Tspect,Fspect,abs(S));
    set(Axes,'YScale','log')
    view(Axes,[0 90]);
    xlabel(Axes,'Time (s)');
    ylabel(Axes,'Frequency (Hz)');
    xlim(Axes,[0 T]);
    ylim(Axes,[fbandlo/3 ftouch*2])
    %titlestr = sprintf('Waveform spectrogram (%.1f,%.1f)',M1,M2);
    %title(titlestr);
    %fname_png = sprintf('Chirp_spectrogram_%d_%d.png',floor(M1),floor(M2));
    %print('-dpng',fname_png);
    return;
    
  end

  % Write audio file after normalizing and tapering to zero at start and after coalescence

  % Prepare audio version of waveform without clipping (max abs value = 1)

  if (RunMode==5)

    haudio = h / max(abs(h)) * 0.9;
    hmax = max(abs(haudio));

  % Define start of audio waveform to play (maximum of 2 seconds in length)
  
    Tplay = 2.0;
    Tmin = max(0,tau-Tplay);
    firstsample = max(ceil(Tmin*fsamp),1);

    % Apply sigmoid taper to start-up of played waveform

    halfperiod_start = min(0.2*Tplay,0.2*tau);
    twindow_start = transpose([0:1/fsamp:halfperiod_start]);
    envelope_start = 0.5-0.5*cos(pi/halfperiod_start*twindow_start);
    for i = 1:size(envelope_start)
      haudio(firstsample+i-1) = haudio(firstsample+i-1)*envelope_start(i);
    end

    % Apply sigmoid taper to end of waveform

    halfperiod_end = 0.01;
    twindow_end = transpose([0:1/fsamp:halfperiod_end]);
    envelope_end = 0.5 + 0.5*cos(pi/halfperiod_end*twindow_end);
    for i = 1:min(length(twindow_end),N-lastsample);
      haudio(lastsample+i) = envelope_end(i)*abs(hmax)*sin(phi(lastsample+i));
    end

    %fname_mp4 = sprintf('Chirp_audio_%d_%d.mp4',floor(M1),floor(M2));
    %audiowrite(fname_mp4,haudio(t>Tmin),fsamp)
    sound(haudio(t>Tmin),fsamp)

    return;
    
  end
  
  % Plot audio waveform value vs time
  
  if (RunMode==6)
    
    cla(Axes);
    plot(Axes,t(firstsample:end), haudio(firstsample:end));
    xlabel('Time (s)');
    ylabel('h(t) [arb. units]');
    %titlestr = sprintf('Audio inspiral waveform in time domain (%.1f,%.1f)',M1,M2);
    %title(titlestr);
    xlim([Tmin T]);
    %fname_png = sprintf('Audio_chirp_waveform_%d_%d.png',floor(M1),floor(M2));
    %print('-dpng',fname_png);
    
    return;
    
  end

  % Create animation of two stars orbiting common center of mass

  % (Define cm as origin and center of animation figure)
  % (Use orbital frequency via Kepler's 3rd law to compute semi-major axis)
  % (The semi-major axis a is the relative separation between the star centers,
  %  where the distances r1 and r2 of the stars from the cm add up to a)

  if (RunMode==1)

    cla(Axes);
    a = freq;
    a = (G * (M1+M2)*Msun ./( (freq/2).^2 * 4*pi^2 ) ).^(1/3);
    fprintf('Initial semi-major axis of relative separation = %e m (%f km)\n',a(1),a(1)/1.e3);
    r1 = M2/(M1+M2)*a;
    r2 = M1/(M1+M2)*a;
    fprintf('Initial r1 for M1 = %f km, r2 for M2 = %f km\n',r1(1)/1.e3,r2(1)/1.e3);

    % Start animation with M1 at phi = 0 (x1 = r1, y1 = 0, x2 = -r2, y2 = 0)
    % At each increment in time, move stars by angle governed by mean frequency
    % during increment but recalculate radii according to frequency
    % at end of increment

    x1 = r1; y1 = r1; x2 = r2; y2 = r2; phi1 = r1; phi2 = r2;
    phi1(1) = 0; phi2(1) = pi;
    for sample = 2:lastsample
      deltaphi = 2*pi * (freq(sample-1)+freq(sample))/2 * dt;
      phi1(sample) = phi1(sample-1) + deltaphi;
      phi2(sample) = phi2(sample-1) + deltaphi;
    end
    x1 = r1.*cos(phi1);
    y1 = r1.*sin(phi1);
    x2 = r2.*cos(phi2);
    y2 = r2.*sin(phi2);

    downsamplefactor = 100;

    % drawnow animation (use a curved "rectangle" as a circle)

    D1 = 2*R1;
    D2 = 2*R2;
    set(Axes,'XScale','linear')
    set(Axes,'YScale','linear')
    circle1 = rectangle(Axes,'Position',[x1(1)-R1 y1(1)-R1 D1 D1],'curvature',[1 1],'FaceColor','red','EdgeColor','red');
    circle2 = rectangle(Axes,'Position',[x2(1)-R2 y2(1)-R2 D2 D2],'curvature',[1 1],'FaceColor','blue','EdgeColor','blue');
    xlim(Axes,[-a(1) a(1)]);
    ylim(Axes,[-a(1) a(1)]);
    xlabel(Axes,' ');
    ylabel(Axes,' ');
    set(Axes,'XTick',[ ]);
    set(Axes,'YTick',[ ]);
    for sample = 2:downsamplefactor:lastsample;
      circle1.Position = [x1(sample)-R1 y1(sample)-R1 D1 D1];
      circle2.Position = [x2(sample)-R2 y2(sample)-R2 D2 D2];
      pause(1.e-6); 
      % drawnow limitrate
    end
    circle1.Position = [x1(lastsample)-R1 y1(lastsample)-R1 D1 D1];
    circle2.Position = [x2(lastsample)-R2 y2(lastsample)-R2 D2 D2];
    return;
    
  end
  
end
