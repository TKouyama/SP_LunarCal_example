;;
;; 観測時間を指定して,その時に見えるはずの月画像をSP月反射率モデルに基づき作成する
;;

;;
;; Required IDL codes
;;
@channel_wavelength.pro
@get_radiance_lism.pro
;@get_geometry_with_spice.pro
@band_response.pro
@lunar_map_plot.pro
@channel_contrib.pro
@f_get_radiance_integ.pro
@band_response.pro

;; センサーの仕様やSPICEカーネルなどを記述したファイルを用意しておく
;; <使い方>
;; read_sp_model_bilinear,,obs_geo, out_irad, out_wav, ofldname = ofldname
;;


pro read_sp_model_bilinear_for_pub $
    ,obs_geo, out_irad, out_wav, out_hyper_image = out_hyper_image, ofldname = ofldname

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
  
  ;;
  ;;
  ;;
  xs = 800; 512*1.5
  ys = 800; 512*1.5
  window_Set,3,xs=xs,ys=ys

  ;; 緯経線を重ねるなら1, 重ねないなら 0 ;;
  ;overplot_latlon = 1
  overplot_latlon = 0

  ;; work_dirにいるようにしておく ;;
  ;work_dir = 'C:\work\ENVI_IDL\LunarCal_pub\'
  ;cd,work_dir

  ;; workディレクトリ名を取得 ;;
  cd,'.',current=current_dir

  ;; 各種必要なディレクトリ名＋ファイル名 ;;
  ;; Dataを収めたディレクトリ ;;
  datadir ='./Parameters/'
  ;; Sol radiacne data ;;
  solfname = datadir+'Gueymard.txt'

  ;; Output ディレクトリ ;;
  if n_elements(ofldname) eq 0 then begin
    ;ofldname =base_dir+'Products\'
    ofldname ='./Outputs/'
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 観測機器の指定, 仕様 ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;
  inst_name = 'ASTER'
  ;inst_name = 'ASTER_BACK'
  ;; カメラ諸パラメータ
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


  ;stop
  ;; obs_geo :: 観測時ジオメトリ情報を収めた構造体。

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
    print,pix_number
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

  wset,0
  if sol_rad_mode eq 'OBS' then begin
    ;; Based on an observation on 2004/04/24 ;;
    ;sol_rad=channel_radiance()/SL_distance_au^2.*1000./!dpi ; W/m2/um/str
    ;plot,wav,sol_rad,title="Solar radiance"
  endif else if sol_rad_mode eq 'LISM' then begin
    ;; Based on Lism standard, Gueymard 2004;;
    get_radiance_lism,sol_rad,solfname=solfname,/test
    sol_rad = sol_rad/SL_distance_au^2.*1000./!dpi
    print,solfname
  endif


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

  loadct,0
  ;window,0,xs=1200,ys=480
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
  ;; new coeff @ 7 ;; G2
  ;k0 = 2.0580020d
  ;k1 = -0.0024387662d
  ;k2 = 1.6992882e-006
  ;k3 = -4.0356982e-010

  ;; new coeff @ 9.57 ;; Hodoyoshi-1
  ;k0 = 2.0884249;     0.014450085
  ;k1 = -0.0025528045;  4.6397498e-005
  ;k2 = 1.8075405e-006;  4.7069029e-008
  ;k3 = -4.3653910e-010;  1.4924942e-011

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

  ;; new coeff @ 31.15 ;; AMICA/NIRS
  ;k0 = 2.0439400
  ;k1 = -0.0023701578
  ;k2 = 1.6708943e-006
  ;k3 = -4.0612406e-010

  ;; new coeff @ 42.74266d ;; OREx
  ;k0 = 2.0859934
  ;k1 = -0.0024645981
  ;k2 = 1.7295874e-006
  ;k3 = -4.1628950e-010

  ;; new coeff @ 45.792776 ;; Hayabusa2 / NIRS3 // VIIRS
  ;k0 = 2.0402810
  ;k1 = -0.0023652289
  ;k2 = 1.6587560e-006
  ;k3 = -4.0193186e-010

  ;; new coeff @ 59.3 ;; Hayabusa2 / ONC
  ;k0 = 1.9915590
  ;k1 = -0.0022089404
  ;k2 = 1.4917248e-006
  ;k3 = -3.4859872e-010

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

    ;; Histogram ;;
    ;result = histogram(tmp_rad* p_coeff[i],binsiz=1,max=200,min=1)
    ;result_x = dindgen(200)+1.5
    ;max_rad_correct[i] = max(result_x[where(result gt 3)])
    ;plot,result_x,result,psym=10
    ;print,max(result_x[where(result gt 5)])
    ;stop

    ;; 160波長すべて出力する
    ;output = reverse(congrid(ccd_rad[*,*,i],60,60),1)
    ;tvscl,reverse(output,2),i

    ;; 地図座標のまま計算 ;;
    map_sp_rad[*,*,i]=sol_rad[i]*map_sp_ref_out[*,*,i]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; 補正係数を輝度値の段階で適応するなら ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ccd_rad[*,*,i]=sol_rad[i] * (ccd_ref[*,*,i] * p_coeff[i])
  endfor
  ;stop

  ;;;;;;;;;;;;;;;;;;;;;
  ;; Irradiance plot ;;
  ;;;;;;;;;;;;;;;;;;;;;
  scaling_factor = atan(dtheta)^2. ;; solid ange is almost equal iFOV^2. for irradiance conversion
  ;print,"##Scale check: ",dtheta, scaling_factor*n_elements(where(tmp_geo gt -!dpi)),n_elements(where(tmp_geo gt -!dpi))

  wset,0
  plot,wav,total_rad*scaling_factor,xrange=[300,1000],xs=1
  oplot,wav,total_rad_correct*scaling_factor,linestyle=2
  wset,0


  ;wset,3
  ;tvscl,ccd_rad[*,*,60]
  ;wset,0

  ;;;;;;;;;;;;
  ;; output ;;
  ;;;;;;;;;;;;
  ;; 補正前データか補正後データか注意する ;;
  ;ofname_csv = 'C:\work\DATA\SELENE_SP\model_results\sample_correct.csv'
  ;ofname_csv = 'C:\work\DATA\Hayabusa2\OREx_Moon\ocams_ega_cal\sample_correct.csv'
  ofname_csv = ofldname+'sample_correct.csv'

  write_csv,ofname_csv,wav,total_rad*scaling_factor,total_rad_correct*scaling_factor $
    ,max_rad_correct,disk_mean_rad_correct,mean_rad_correct $
    ,header = ['Wavelength','Raw irradiance','Corrected irradiance','Max radiance','Dayside mean radiance','Disk mean radiance']

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 全波長分データを出力する ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;ofname = ofldname+'MODEL_rad_band_full.fits'
  ;print,"Output: ",ofname
  ;writefits,ofname,float(ccd_rad)

  ;;;;;;;;;;;;
  ;; 論文用 ;;
  ;;;;;;;;;;;;
  window,4,xs=720,ys=360
  loadct,0,/silent
  ;; wav[40] = 752.8nm

  tvscl,congrid(map_sp_rad[*,*,40]<100,720,360)
  plots,obs_geo.ssc_longitude*2.+360,obs_geo.ssc_latitude*2.+180,/device,psym=1,symsize=3,color=255,thick=2
  plots,obs_geo.ssl_longitude*2.+360,obs_geo.ssl_latitude*2.+180,/device,psym=7,symsize=3,color=255,thick=2

  ;; reflectance ;;
  ;tvscl,congrid(map_sp_ref_ori[*,*,40]<0.25,720,360)

  ccd_siz = size(ccd_rad)
  ;plot,wav,ccd_rad[ccd_siz[1]/2,ccd_siz[2]/2,*],xtitle='Wavelength (nm)',ytitle='Radiance (W/m2/sr/um)'
  mean_r = dblarr(160)
  for mean_i=0, 160-1, 1 do begin
    tmp_ccd = ccd_rad[*,*,mean_i]
    mean_r[mean_i] = mean(tmp_ccd[where(tmp_ccd gt 0)])
  endfor
  plot,wav,mean_r,xtitle='Wavelength (nm)',ytitle='Radiance (W/m2/sr/um)'
  ;stop

  loadct,39,/silent
  wait,0.001

  wset,0

  ;; --- ここまでで波長ごとの計算終了  --- ;;

  print,"Output...: "
  ;out_irad = total_rad*scaling_factor
  out_irad = total_rad_correct*scaling_factor
  out_wav = wav


  return
end


