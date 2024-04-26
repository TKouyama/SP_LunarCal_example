;;
;; ASTERのカメラスペックを読み込む Backもあり。
;;
function get_cam_info_aster,inst_name=inst_name
  if n_elements(inst_name) eq 0 then begin
    dtheta = 21.3e-6
    pix_number=[640,640] ;; 仮想CCD1辺のピクセル数
  endif else begin
    if inst_name eq 'aster' or inst_name eq 'ASTER' then begin
      ;; Aster Nadir ifov ~ 21.3 mrad/pix (光軸中心における1pixあたりの見込み角)
      dtheta = 21.3e-6 ;; from JAROS
      pix_number=[640,640] ;; 仮想CCD1辺のピクセル数
    endif else if inst_name eq 'aster_back' or inst_name eq 'ASTER_BACK' then begin
      ;; Aster Back ifov ~ 18.6 mrad/pix (光軸中心における1pixあたりの見込み角)
      dtheta = 18.6e-6 ;; from JAROS
      pix_number=[640,640] ;; 仮想CCD1辺のピクセル数
    endif else if inst_name eq 'human' then begin
      ;; Human ~ 60 arcsec/pix
      dtheta = 60./3600.*!dpi/180. ;; ifov (radian) (i.e. 光軸中心における1pixあたりの見込み角)
      pix_number=[64,64] ;; 仮想CCD1辺のピクセル数
    endif else if inst_name eq 'test' then begin
      dtheta = 2./3600.*!dpi/180. ;; ifov (radian) (i.e. 光軸中心における1pixあたりの見込み角)
      pix_number=[512,512] ;; 仮想CCD1辺のピクセル数
    endif else begin
      ;; 仮のセンサとしてASTERを与える ;;
      dtheta = 21.3e-6
      pix_number=[512,512] ;; 仮想CCD1辺のピクセル数
    endelse
    
  endelse

  cam_info = dblarr(3)
  cam_info[0] = dtheta
  cam_info[1] = pix_number[0]
  cam_info[2] = pix_number[1]

  return,cam_info
end