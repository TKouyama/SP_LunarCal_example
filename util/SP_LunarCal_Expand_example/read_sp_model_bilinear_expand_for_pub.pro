;;
;; Core program for simulating a Moon observation
;; In this code, spatial resolution is set to be same as that for ASTER (15m GSD from 700 km)
;; t.kouyama@aist.go.jp
;;

;;
;; <usage>
;; read_sp_model_bilinear,obs_geo $
;;                      , out_wav ,out_hyper_image, out_irad
;;

;;
;; Input: obs_geo: a strucutre for containing observation geometry
;;


pro read_sp_model_bilinear_expand_for_pub $
    ,obs_geo, out_wav, out_hyper_image, out_irad, out_ccd_geo, datadir = datadir

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Geometry情報を更新したいとき ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if n_elements(obs_geo) eq 0 then begin
    AU = 149597870700d0/1000d0 ;; 1 AU (km)
    ;; ターゲットまでの距離 ;;
    mean_MD = 384400d

    SL_distance = AU
    LT_distance = mean_MD
    ;; for output ;;
    NA_deg = -90;; General case

    sub_solar_lon_deg = -5.98037
    sub_solar_lat_deg = 0d
    sub_terra_lon_deg = -12.0
    sub_terra_lat_deg = -8.0
    phase_deg = abs(-10d)

    ;; 今回は使わないパラメータ ;;
    phi_deg = 0d
    lam_deg = 0d

    act_Rm = 3474.3/2. ;; km
    ang_diam_deg = 2.*asin(act_Rm/LT_distance)*180./!dpi

    tags = ['ST_distance','TO_distance','ssl_longitude','ssl_latitude' $
      ,'ssc_longitude', 'ssc_latitude', 'N_azimuth', 'phi', 'lam','phase']
    obs_geo = create_struct(tags $
      ,SL_distance $
      ,LT_distance $
      ,sub_solar_lon_deg $
      ,sub_solar_lat_deg $
      ,sub_terra_lon_deg $
      ,sub_terra_lat_deg $
      ,NA_deg $
      ,phi_deg $
      ,lam_deg $
      ,phase_deg $
      )    
  endif
  

  ;; 緯経線を重ねるなら1, 重ねないなら 0 ;;
  ;overplot_latlon = 1
  overplot_latlon = 0

  ;; currentディレクトリ名を取得 ;;
  cd,'./',current=current_dir

  ;; 各種必要なディレクトリ名＋ファイル名 ;;
  ;; Dataを収めたディレクトリ ;;
  if n_elements(datadir) eq 0 then begin
    datadir ='./Parameters/'    
  endif
  ;; Sol radiacne data ;;
  solfname = datadir+'Gueymard.txt'

  ;; Output ディレクトリ ;;
  ;if n_elements(ofldname) eq 0 then begin
  ;  ofldname ='./Outputs/'
  ;endif

  ;;;;;;;;;;;;;;;;;;;;;;
  ;; 観測機器の指定, 仕様 ;;
  ;;;;;;;;;;;;;;;;;;;;;;
  ;; カメラ諸パラメータ
  ;; Simulation image size
  xs = 800; 512*1.5
  ys = 800; 512*1.5

  ;;
  ;; number of bands
  ;;
  ;n_bands = 160l
  n_bands = 204l

  ;; spatial resolution
  inst_name = 'ASTER'
  cam_info = get_cam_info_aster(inst_name=inst_name) ; or (inst_name='ASTER_BACK')

  dtheta = cam_info[0]
  ;pix_number = [cam_info[1],cam_info[2]]
  ;pix_number = [800,800]
  pix_number = [xs,ys]
  ccd_center = (pix_number-1.)/2.

  ;;;;;;;;;;;;;;;;;;;;
  ;; 観測条件を決定 ;;
  ;;;;;;;;;;;;;;;;;;;;
  integ_zoom_ratio = 1. ;; 整数がよい, ASTER条件では2が最大

  print,"!!! Intergrating zoom ratio: ",integ_zoom_ratio
  print,"!!! NA_deg: ",obs_geo.N_azimuth
  print,"!!! IFOV: ",dtheta

  if integ_zoom_ratio gt 1 then begin
    dtheta = dtheta/integ_zoom_ratio
    pix_number=[pix_number[0]*integ_zoom_ratio,pix_number[1]*integ_zoom_ratio]
    ;ccd_center = ccd_center*integ_zoom_ratio
    ccd_center = (pix_number-1.)/2.

  endif

  ;;
  ;; 画像サイズが大きすぎるときは止める ;;
  ;;
  if pix_number[0] gt 2048 then begin
    print,"Check image frame size: ", pix_number
    stop
  endif


  ;; Phase angle new ;;
  tmp_obs = [cos(obs_geo.ssc_latitude*!dpi/180.)*cos(obs_geo.ssc_longitude*!dpi/180.) $
    ,cos(obs_geo.ssc_latitude*!dpi/180.)*sin(obs_geo.ssc_longitude*!dpi/180.) $
    ,sin(obs_geo.ssc_latitude*!dpi/180.)]
  tmp_sol = [cos(obs_geo.ssl_latitude*!dpi/180.)*cos(obs_geo.ssl_longitude*!dpi/180.) $
    ,cos(obs_geo.ssl_latitude*!dpi/180.)*sin(obs_geo.ssl_longitude*!dpi/180.) $
    ,sin(obs_geo.ssl_latitude*!dpi/180.)]
  tmp_phase_angle = acos(total(tmp_obs*tmp_sol))*180./!dpi
  print,"Phase angle tmp: ",tmp_phase_angle


  ;help,obs_geo ;; ← 中身を確認したい時
  ;stop

  ;; 太陽までの距離 ;;
  AU = 149597870700d0/1000d0 ;; 1 AU (km)
  SL_distance_au = obs_geo.ST_distance/AU

  ;; ターゲットまでの距離 ;;
  mean_MD = 384400d
  ;obs_geo.TO_distance = mean_MD
  LT_distance = obs_geo.TO_distance

  act_Rm = 3474.3/2. ;; km
  ang_diam_deg = 2.*asin(act_Rm/LT_distance)*180./!dpi

  print,"Distance: ", SL_distance_au
  print,"Distance: ",LT_distance, " [km]"
  print,"Ang_diam: ", ang_diam_deg, " [deg]"
  print,"Integ zoom", integ_zoom_ratio
  ;stop

  ;; spモデル波長の取得
  wav = channel_wavelength()
  ;print,wav[40] ; <= 752.8 nm

  ;; 入射光の強さ
  ;; フラックスなので距離補正が必要
  sol_rad_mode = 'LISM'

  ;wset,0
  ;; Based on Lism standard, Gueymard 2004;;
  get_radiance_lism,sol_rad,solfname=solfname;,/test
  sol_rad = sol_rad/SL_distance_au^2.*1000./!dpi
  print,solfname

  ;;;;;;;;;;;;;;;;;;;;
  ;; モデル画像作成 ;;
  ;;;;;;;;;;;;;;;;;;;;
  ;; 2次元アレイCCD撮像を想定したモデル画像を作る ;;
  ;; CCDモデル画像では(pix_number-1)/2を投影の中心点とする(ex. pix_number = 512の時は255.5)
  ;ccd_center = (pix_number-1.)/2.

  ;; 波長ごとのRefrectanceを読み込んでccdに写る形に投影
  ;; 投影条件そのままで出力するので出力配列は上下左右反転している。必要に応じてさらに反転させる
  ;; Bilinear interpolationを利用している
  ;; get_geometry_with_spice で取得するobs_geo構造体が必要

  ;; Reddign 対策なし ;;
  lunar_map_plot,obs_geo,ccd_ref $
    ,ccd_geo=ccd_geo,dtheta=dtheta,pix_number=pix_number,ccd_center=ccd_center $
    ,datadir=datadir,correct_reflectance=correct_reflectance $
    ,map_sp_ref_out=map_sp_ref_out,map_sp_ref_ori=map_sp_ref_ori

  ;; 反射率に太陽光輝度値をかけて反射光輝度値に変換
  cube_siz = size(ccd_ref)
  wav_number = cube_siz[3]

  total_rad = dblarr(wav_number)
  total_rad_correct = dblarr(wav_number)
  mean_rad_correct = dblarr(wav_number)
  disk_mean_rad_correct = dblarr(wav_number)
  max_rad_correct = dblarr(wav_number)

  ;; lat ;;
  tmp_geo = ccd_geo[*,*,1]
  tmp_geo_pos = where(tmp_geo gt -!dpi,geo_count)

  ;; emission ;;
  tmp_emi = ccd_geo[*,*,2]
  tmp_emi_pos = where(tmp_emi lt 75.*!dpi/180.)

  ;; 補正パラメータ初期化
  k0 = 1d
  k1 = 0d
  k2 = 0d
  k3 = 0d

  ;; ROLO correct coefficient
  ;; in paper @ 27.7 ;;
  ;k0 = 1.885
  ;k1 = -1.891e-3
  ;k2 = 1.257e-6
  ;k3 = -2.984e-10
  ;tmp_p_coeff0 = k0+k1*wav+k2*wav^2.+k3*wav^3.

  ;; new coeff @ -27.7 ;;
  k0 = 2.1453804
  k1 = -0.0026764001
  k2 = 1.9551764e-006
  k3 = -4.8732791e-010

  p_coeff = k0*replicate(1,n_elements(wav))+k1*wav+k2*wav^2.+k3*wav^3.

  ;;;;;;;;;;;;;;;;
  ;; 適応しない ;; 月画像を介した比較の時はこちらが良い?
  ;;;;;;;;;;;;;;;;
  ;p_coeff[*] = 1.

  ;; 初期化 ;;
  ccd_rad = dblarr(pix_number[0],pix_number[1],wav_number)
  map_sp_rad = map_sp_ref_out*0
  for i=0, wav_number-1,1 do begin
    ;; 反射率をRadianceに変換する ;;
    ccd_rad[*,*,i]=sol_rad[i]*ccd_ref[*,*,i]

    ;; Evaluation of Irradiance ;;
    tmp_rad = ccd_rad[*,*,i]
    tmp_total_rad = total(tmp_rad[where(tmp_rad gt 0)])

    ;; 補正前データ ;;
    total_rad[i] = tmp_total_rad
    ;; 補正後データ ;;
    total_rad_correct[i] = tmp_total_rad * p_coeff[i]

    ;; 影も含めた平均 ;;
    mean_rad_correct[i] = total_rad_correct[i]/geo_count

    ;; 日照面だけの平均 ;;
    disk_mean_rad_correct[i] = mean(tmp_rad[where(tmp_rad gt 0)])* p_coeff[i]
    ;print,mean_rad_correct[i], mean(tmp_rad[where(tmp_rad gt 0)])* p_coeff[i]

    ;; 輝度最大値 ;;
    max_rad_correct[i] = max(tmp_rad[tmp_emi_pos]*p_coeff[i])

    ;; 地図座標のまま計算 ;;
    map_sp_rad[*,*,i]=sol_rad[i]*map_sp_ref_out[*,*,i]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; 補正係数を輝度値の段階で適応するなら ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ccd_rad[*,*,i]=sol_rad[i] * (ccd_ref[*,*,i] * p_coeff[i])
  endfor
  ;stop

  ;;
  ;; Irradiance ;;
  ;;
  scaling_factor = atan(dtheta)^2. ;; solid ange is almost equal iFOV^2. for irradiance conversion
  ;print,"##Scale check: ",dtheta, scaling_factor*n_elements(where(tmp_geo gt -!dpi)),n_elements(where(tmp_geo gt -!dpi))

  ;;
  ;; Mean radiance
  ;;
  ccd_siz = size(ccd_rad)
  mean_r = dblarr(n_bands)
  for mean_i=0, n_bands-1, 1 do begin
    tmp_ccd = ccd_rad[*,*,mean_i]
    mean_r[mean_i] = mean(tmp_ccd[where(tmp_ccd gt 0)])
  endfor
  ;plot,wav,ccd_rad[ccd_siz[1]/2,ccd_siz[2]/2,*],xtitle='Wavelength (nm)',ytitle='Radiance (W/m2/sr/um)'

  ;;
  ;; QL Examples
  ;;
  ;wset,0
  ;plot,wav,total_rad*scaling_factor,xrange=[300,1000],xs=1
  ;oplot,wav,total_rad_correct*scaling_factor,linestyle=2
  ;wset,0

  ;wset,3
  ;tvscl,ccd_rad[*,*,60]
  ;wset,0

  ;loadct,0
  ;window,0,xs=1200,ys=480
  ;window,4,xs=720,ys=360
  ;loadct,0,/silent
  ;; simulated image at a certain band ;; ;; wav[40] = 752.8nm
  ;tvscl,congrid(map_sp_rad[*,*,40]<100,720,360)
  ;plots,obs_geo.ssc_longitude*2.+360,obs_geo.ssc_latitude*2.+180,/device,psym=1,symsize=3,color=255,thick=2
  ;plots,obs_geo.ssl_longitude*2.+360,obs_geo.ssl_latitude*2.+180,/device,psym=7,symsize=3,color=255,thick=2

  ;; reflectance map ;;
  ;tvscl,congrid(map_sp_ref_ori[*,*,40]<0.25,720,360)

  ;; mean radiance ;;
  ;plot,wav,mean_r,xtitle='Wavelength (nm)',ytitle='Radiance (W/m2/sr/um)'

  ;loadct,39,/silent
  ;wait,0.001
  ;wset,0

  ;; --- ここまでで波長ごとの計算終了  --- ;;

  print,"Done...: "
  ;out_irad = total_rad*scaling_factor
  ;; integrated irradiance ;;
  out_irad = total_rad_correct*scaling_factor
  ;; channel wavelength
  out_wav = wav
  ;; simulated 2D image (hyperspectral)
  out_hyper_image = ccd_rad
  
  out_ccd_geo = ccd_geo

  return
end


