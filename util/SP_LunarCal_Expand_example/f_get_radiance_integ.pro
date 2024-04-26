;;
;; obs_wavからsel_wavへ変更するために(台形)積分するプログラム
;; 波長刻みが等間隔でないときは不正確の可能性(2022.01.16)
;;
function f_get_radiance_integ,obs_wav, obs_rad,sel_wav
  
  obs_wav_p = shift(obs_wav,-1)
  obs_wav_m = shift(obs_wav,1)
  obs_int_p = (obs_wav_p-obs_wav)
  obs_int_m = (obs_wav-obs_wav_m)
  obs_int = obs_int_p + obs_int_m
  n_data = n_elements(obs_wav)
  
  siz = size(sel_wav)
  
  ;; 一つサイズの大きい波長配列を用意
  tmp_sel_wav = dblarr(siz[1]+1)
  tmp_sel_wav[0:siz[1]-1]=sel_wav
  tmp_sel_wav[siz[1]] = (sel_wav[siz[1]-1]-sel_wav[siz[1]-2])+sel_wav[siz[1]-1]
  
  sel_rad = dblarr(siz[1])
  position = 0l
  position1 = 0l
  
  ;; wav 0------wav 1------wav 2
  ;;        |<--range-->|
  for i=0, siz[1]-1,1 do begin
    if i eq 0 then begin
      while(obs_wav[position] le tmp_sel_wav[i]) do begin
        if position ge n_data-1 then break
        position++
      endwhile
      
      position1 = position
      while(obs_wav[position1] le (tmp_sel_wav[i]+tmp_sel_wav[i+1])/2.) do begin
        if position1 ge n_data-1 then break
        position1++
      endwhile
    endif else begin
      while(obs_wav[position] le (tmp_sel_wav[i]+tmp_sel_wav[i-1])/2.) do begin
        if position ge n_data-1 then break
        position++
      endwhile
      
      position1 = position
      while(obs_wav[position1] le (tmp_sel_wav[i]+tmp_sel_wav[i+1])/2.) do begin
        if position1 ge n_data-1 then break
        position1++
      endwhile
    endelse
    
    if position ge n_data then return, sel_rad
    ;; 単純平均 ;;
    ;sel_rad[i]=mean(obs_rad[position:position1])
    
    ;; 台形積分平均 ;;
    if position1 eq position then begin
      if position1 eq 0 then begin
        sel_rad[i]=0d
      endif else if position eq n_data-1 then begin
       sel_rad[i]=0d
      endif else begin
       sel_rad[i]=mean(obs_rad[position:position1])
      endelse
    endif else if position1 gt position then begin
      sel_rad[i]=0d
      tmp_total_wav = 0d

      for j = position, position1-1, 1 do begin
        d_wav = (obs_wav[j+1]-obs_wav[j])/1000d
        sel_rad[i]+=(obs_rad[j+1]+obs_rad[j])*d_wav/2d
        tmp_total_wav += d_wav
      endfor
      if tmp_total_wav eq 0 then begin
        print, i,position,position1
        stop
      endif
      sel_rad[i]/= tmp_total_wav

    endif

    ;; Linear interpolation
    ;min_pos = position-1
    ;a = obs_wav[min_pos+1]-sel_wav[i]
    ;b = sel_wav[i]-obs_wav[min_pos]
    ;sel_rad[i] = (a*obs_rad[min_pos+1] + b*obs_rad[min_pos])/(a+b)
  endfor
  
  return,sel_rad
end
