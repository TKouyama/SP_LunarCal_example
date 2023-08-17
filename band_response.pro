;;
;; Aster Band Response
;;

function interpolate_band_response,sb_wl,sb_re,ch_wav,int_re,contrib
  ;ch_wav = channel_wavelength(arr_size=arr_size)/1000.

  arr_size=1000
  ch_wav = dindgen(arr_size)/arr_size*(1000.-400.)/1000.+0.4

  ch_wav_m = shift(ch_wav,1)
  ch_wav_p = shift(ch_wav,-1)
  d_ch_wav = (ch_wav_p-ch_wav_m)/2.
  d_ch_wav[0]=0
  d_ch_wav[arr_size-1]=0

  ;ch_wav_pos = where(ch_wav gt min(sb_wl-0.1) and ch_wav lt max(sb_wl+0.1))
  ch_wav_pos = where(ch_wav gt min(sb_wl-0.01) and ch_wav lt max(sb_wl+0.01))
  ch_wav = ch_wav[ch_wav_pos]
  d_ch_wav = d_ch_wav[ch_wav_pos]
  int_re = INTERPOL( sb_re, sb_wl, ch_wav)>0

  ;; 寄与率に変化する ;;
  ;total_S = total(d_ch_wav)
  ;contrib = int_re*d_ch_wav/total_S
  tmp_contrib = int_re*d_ch_wav
  total_S = total(tmp_contrib)
  
  ;; Test ;;
  contrib = tmp_contrib/total_S
  
  contrib = int_re
  c_arr_size = size(contrib)
  return,c_arr_size[1]
end

;;
;;
;;
pro band_response,ch_wav_all=ch_wav_all, contrib_all = contrib_all,no_plot=no_plot
  ifldname = 'C:\work\DATA\ASTER_LunaCal\BandResponse\'
  ifname = 'vnir.txt'
  ofname = 'vnir_tra.txt'

  cd, ifldname, current=old
  
  line=''
  count=0
  openr,1,ifname
  ;; 4行読み飛ばし
  readf,1,line
  readf,1,line
  readf,1,line
  readf,1,line

  ;; 作業用変数
  tmp_band1_wl = 0.
  tmp_band2_wl = 0.
  tmp_band3N_wl = 0.
  tmp_band3B_wl = 0.
  
  tmp_band1_re = 0.
  tmp_band2_re = 0.
  tmp_band3N_re = 0.
  tmp_band3B_re = 0.
  
  band1_wl=dblarr(100)
  band2_wl=dblarr(100)
  band3N_wl=dblarr(100)
  band3B_wl=dblarr(100)

  band1_re=dblarr(100)
  band2_re=dblarr(100)
  band3N_re=dblarr(100)
  band3B_re=dblarr(100)

  while(not eof(1) and count lt 100) do begin
    readf,1,tmp_band1_wl,tmp_band1_re $
           ,tmp_band2_wl,tmp_band2_re $
           ,tmp_band3N_wl,tmp_band3N_re $
           ,tmp_band3B_wl,tmp_band3B_re

    band1_wl[count]=tmp_band1_wl
    band2_wl[count]=tmp_band2_wl
    band3N_wl[count]=tmp_band3N_wl
    band3B_wl[count]=tmp_band3B_wl

    band1_re[count]=tmp_band1_re
    band2_re[count]=tmp_band2_re
    band3N_re[count]=tmp_band3N_re
    band3B_re[count]=tmp_band3B_re

    count++
  endwhile
  close,1

  band1_wl = band1_wl[0:count-1]
  band2_wl = band2_wl[0:count-1]
  band3N_wl = band3N_wl[0:count-1]
  band3B_wl = band3B_wl[0:count-1]

  band1_re = band1_re[0:count-1]
  band2_re = band2_re[0:count-1]
  band3N_re = band3N_re[0:count-1]
  band3B_re = band3B_re[0:count-1]

  sb_wl = band1_wl ;; Selected band wave length
  sb_re = band1_re ;; Selected band response
  siz1 = interpolate_band_response(sb_wl,sb_re,ch_wav1,int_re1,contrib1)

  sb_wl = band2_wl ;; Selected band wave length
  sb_re = band2_re ;; Selected band response
  siz2 = interpolate_band_response(sb_wl,sb_re,ch_wav2,int_re2,contrib2)
  
  sb_wl = band3N_wl ;; Selected band wave length
  sb_re = band3N_re ;; Selected band response
  siz3N = interpolate_band_response(sb_wl,sb_re,ch_wav3N,int_re3N,contrib3N)

  ch_wav_all = [ch_wav1,ch_wav2,ch_wav3N]
  contrib_all = [contrib1,contrib2,contrib3N]
  
  ;; Band 3B ;;
  sb_wl = band3B_wl ;; Selected band wave length
  sb_re = band3B_re ;; Selected band response
  siz3B = interpolate_band_response(sb_wl,sb_re,ch_wav3B,int_re3B,contrib3B)

  ;ch_wav_all = [ch_wav_all,ch_wav3B]
  ;contrib_all = [contrib_all,contrib3B]

  if keyword_set(no_plot) then begin
  endif else begin
    plot,ch_wav_all,contrib_all,/nodata
    oplot,ch_wav1,contrib1,linestyle=0
    oplot,ch_wav2,contrib2,linestyle=3
    oplot,ch_wav3N,contrib3N,linestyle=2
    oplot,ch_wav3B,contrib3B,linestyle=1
    print,total(contrib3N)

    ofname = 'vnir_band1.csv'
    close,11
    openw,11,ofname
    for i=0, siz1-1,1 do begin
      printf,11,ch_wav1[i],",",contrib1[i]
    endfor
    close,11
  endelse


  cd, old
end