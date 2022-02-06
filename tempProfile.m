
% Ez nagyjabol azonos feltetelekkel keszult, mint a Wakid Tplate=48, Tamb=20 kiserlet
% Tort drapp vonal: Wakid.
% Cakkos kek: Croft modell a pontosabbnak mondott szamitassal az also tartomanyban - a ket esetet nyilvan nem hozta ossze.
% Folytonos narancs: Croft a kozelito szamitassal.
% Ez utobbi jo kozelites lesz, mert mire elter a Wakid-fele mert ertekektol, mar alig van fenytores a kicsi homerseklet-gradiens miatt.
% A gorbe elejerol szerintem nem sok mondhato a Wakid-kiserlet keves pontja miatt.
function result = tempProfile(aHeight, aExact) % meters
  Tplate = 326.7;            % Kelvin
  heightLimit = 0.25;        % cm
  B = 0.041;
  delta = 1.4;
  Tamb = 297.6;
  Theta0 = 0.03;
  Z0 = 0.129;
  Z = aHeight * 100;
  if(Z > heightLimit)
    result = Tamb * exp(B * Z .^ (1 - delta));
  else
    if(aExact)
      result = Tamb * exp(Theta0 * 3.1 * sqrt(Z0 / Z));
    else
      result = Tamb * exp(B * Z .^ -0.2);
    endif
  endif
endfunction
