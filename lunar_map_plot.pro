;;
;; 衛星直下点緯度経度・太陽直下点緯度経度を与えて
;; 観測を再現するプログラム.
;; 各画素の輝度値は月面代表点によって算出。(面積分計算を適応するか否か)
;; 衛星による観測は有限距離からであることを考慮、太陽は無限遠近似.
;; 月は取り急ぎ球体近似.(高度や傾斜未考慮)
;; 光時・光路補正は直下点緯度経度を求める際に行っているものとする.
;;
;; Albego groupは0.5x0.5°刻みなのであらかじめマップで反射率換算し、Bilinearで補間する.
;; 解像度が低い観測に対応するには高解像度画像を作ってビニング処理する
;;
;; 参考文献
;; Lunar lamert lawはYokota et al., 2011, Icarus
;;

;; 解全要素チェック ;;
;; lunar_map_plot 内 ;;
;; ccd_geo をうまく更新してやれば歪み考慮後の画像に適応できる? ;; 2014/04/04 ;;
;; mapの切れ目対策を追記 ;; 2015/10/07
;;


Forward_Function latlon_map
Forward_Function map_ccd_latlon
Forward_Function map_read_photometric_params
Forward_Function map_photomertic_function
Forward_Function read_sp_reflectance_rcc

;; 地図投影上の位相角関係
function latlon_map,lat,lon,crange,Rv,f_Lv,f_Nv,ssl_lat,ssl_lon,ssc_lat,ssc_lon $
                   ,c_map_distance,c_map_inc,c_map_emi,c_map_phase

    ;; 太陽直下点ベクトル
    l_ssl = dblarr(3)
    l_ssl[0] = cos(ssl_lat)*cos(ssl_lon)
    l_ssl[1] = cos(ssl_lat)*sin(ssl_lon)
    l_ssl[2] = sin(ssl_lat)

    ;; 衛星直下点
    l_ssc = dblarr(3)
    l_ssc[0] = cos(ssc_lat)*cos(ssc_lon)
    l_ssc[1] = cos(ssc_lat)*sin(ssc_lon)
    l_ssc[2] = sin(ssc_lat)

    R2 = Rv^2.
    D2 = Crange^2.
    DR2 = 2.*Rv*Crange

    ;; 各緯度経度の規格化された位置
    l_map_x = cos(lat)*cos(lon)
    l_map_y = cos(lat)*sin(lon)
    l_map_z = sin(lat)

    ;; Incident Angle
    tmp_prod = l_ssl[0]*l_map_x+l_ssl[1]*l_map_y+l_ssl[2]*l_map_z
    c_map_inc = acos(tmp_prod *(abs(tmp_prod) le 1.))*(abs(tmp_prod) le 1.)

    ;; Emission Angle ;; 衛星直下点で定義域を出すので対策した (2012/06/15)
    tmp_prod = l_ssc[0]*l_map_x+l_ssc[1]*l_map_y+l_ssc[2]*l_map_z
    tmp_angle = acos((tmp_prod<1)>(-1))
    ;tmp_angle = acos(tmp_prod *(abs(tmp_prod) le 1.))*(abs(tmp_prod) le 1.)
    tmp_X = sqrt((R2+D2-DR2*cos(tmp_angle))>0.)
    tmp_angle2 = acos( (((R2+tmp_X^2.-D2)/(2.*Rv*tmp_X))<1)>(-1) )
    c_map_emi = !dpi-tmp_angle2

    ;; Bright
    c_map_imu = cos(c_map_inc) > 0.
    c_map_emu = cos(c_map_emi) > 0.

    ;; Phase angle
    ;; 衛星からみた天体基本3軸
    f_E1 = crossp(f_Lv, f_Nv)/cos(ssc_lat)
    f_E0 = crossp(f_E1, f_Nv)
    f_E2 = f_Nv

    ;; 衛星から見た天体表面の方向
    V_center = f_Lv*Crange
    map_Lv_x = Rv*(l_map_x*f_E0[0] + l_map_y*f_E1[0] + l_map_z*f_E2[0]) + V_center[0]
    map_Lv_y = Rv*(l_map_x*f_E0[1] + l_map_y*f_E1[1] + l_map_z*f_E2[1]) + V_center[1]
    map_Lv_z = Rv*(l_map_x*f_E0[2] + l_map_y*f_E1[2] + l_map_z*f_E2[2]) + V_center[2]
    c_map_distance = sqrt(map_Lv_x^2.+map_Lv_y^2.+map_Lv_z^2.)
    map_Lv_x/=c_map_distance
    map_Lv_y/=c_map_distance
    map_Lv_z/=c_map_distance

    ;; 行列化
    r_mat = dblarr(3,3)
    r_mat[0,*]=f_E0
    r_mat[1,*]=f_E1
    r_mat[2,*]=f_E2
    i_r_mat = invert(r_mat)

    ;; ssl lon - ssc lon ベクトル
    tmp_l_ssl = dblarr(3)
    tmp_l_ssl[0] = cos(ssl_lat)*cos(ssl_lon-ssc_lon)
    tmp_l_ssl[1] = cos(ssl_lat)*sin(ssl_lon-ssc_lon)
    tmp_l_ssl[2] = sin(ssl_lat)
    
    ;tmp_l_ssl[0] = cos(ssl_lat-ssc_lat)*cos(ssl_lon-ssc_lon)
    ;tmp_l_ssl[1] = cos(ssl_lat-ssc_lat)*sin(ssl_lon-ssc_lon)
    ;tmp_l_ssl[2] = sin(ssl_lat-ssc_lat)
    ;; Rotate
    ssl_vec = r_mat##tmp_l_ssl
    ssl_vec /= sqrt(total(ssl_vec^2.))

    ;; Phase angle
    c_map_phase = acos(-map_Lv_x*ssl_vec[0] $
                       -map_Lv_y*ssl_vec[1] $
                       -map_Lv_z*ssl_vec[2])
    return,0  
end

;; CCD各座標にLon Latを対応付ける
;; 引数に ccd_centerを加えた(2012/06/22)
function map_ccd_latlon,dtheta,pix_number,ccd_center,crange,Rv,f_Lv,f_Nv,ssl_lat,ssl_lon,ssc_lat,ssc_lon

    ;; CCD 座標定義
    c_ccd_x = dblarr(pix_number[0],pix_number[1])
    for j=0,pix_number[1]-1,1 do begin
      c_ccd_x[*,j] = dindgen(pix_number[0])
    endfor
    c_ccd_y = dblarr(pix_number[0],pix_number[1])
    for i=0,pix_number[0]-1,1 do begin
      c_ccd_y[i,*] = dindgen(pix_number[1])
    endfor
    
    ;;;;;;;;;;;;
    ;; CCD中心 ;;
    ;;;;;;;;;;;;
    ccd_center_x = ccd_center[0]
    ccd_center_y = ccd_center[1]    
    
    ;; 各CCD素子の視線ベクトル
    c_ccd_Lv = dblarr(pix_number[0],pix_number[1],3)
    c_ccd_Lv[*,*,0] = -(c_ccd_x-ccd_center_x)*dtheta
    c_ccd_Lv[*,*,1] = -(c_ccd_y-ccd_center_y)*dtheta
    c_ccd_Lv[*,*,2] = 1.
    tmp_c_ccd_Lv = c_ccd_Lv
    for i=0,2,1 do begin
      tmp_c_ccd_Lv[*,*,i] /= sqrt(c_ccd_Lv[*,*,0]^2.+c_ccd_Lv[*,*,1]^2.+c_ccd_Lv[*,*,2]^2.)
    endfor
    c_ccd_Lv = tmp_c_ccd_Lv

    ;; Lvベクトル
    if n_elements(f_Lv) eq 0 then f_Lv = [0,0,1]
    
    ;; 各CCDから出る視線ベクトルとLvベクトルのなす角
    c_ccd_Lv_ang = acos(c_ccd_Lv[*,*,0]*f_Lv[0]+c_ccd_Lv[*,*,1]*f_Lv[1]+c_ccd_Lv[*,*,2]*f_Lv[2])
    ;; リムとLvベクトルのなす角 (= 最大角度)
    limit_ang  = asin(Rv/Crange)


    ;; 天体表面の投影点までの距離
    sq_surface_dis = Rv^2.-Crange^2.*sin(c_ccd_Lv_ang)^2.
    surface_dis = (Crange*cos(c_ccd_Lv_ang)-sqrt(sq_surface_dis > 0))*(c_ccd_Lv_ang le Limit_ang)

    ;; 衛星から見たときの天体表面座標-天体中心座標=天体中心から見た表面座標(ただし衛星座標系の基底表現)
    V_center = f_Lv*Crange
    V_sx = (c_ccd_Lv[*,*,0]*surface_dis - V_center[0])*(c_ccd_Lv_ang le Limit_ang)
    V_sy = (c_ccd_Lv[*,*,1]*surface_dis - V_center[1])*(c_ccd_Lv_ang le Limit_ang)
    V_sz = (c_ccd_Lv[*,*,2]*surface_dis - V_center[2])*(c_ccd_Lv_ang le Limit_ang)
    abs_V_s = sqrt(V_sx^2.+V_sy^2.+V_sz^2.)

    pos = where(abs_V_s ne 0)
    ;print,mean(abs_V_s[pos])

    V_sx[pos]/=abs_V_s[pos]
    V_sy[pos]/=abs_V_s[pos]
    V_sz[pos]/=abs_V_s[pos]

    ;; 基本3軸
    f_E1 = crossp(f_Lv, f_Nv)/cos(ssc_lat)
    f_E0 = crossp(f_E1, f_Nv)
    f_E2 = f_Nv
    
    ;; 行列化
    r_mat = dblarr(3,3)
    r_mat[0,*]=f_E0
    r_mat[1,*]=f_E1
    r_mat[2,*]=f_E2
    i_r_mat = invert(r_mat)

    ;; 座標→緯度経度
    tmp_vec = dblarr(3)
    c_ccd_lon = dblarr(pix_number[0],pix_number[1])
    c_ccd_lat = dblarr(pix_number[0],pix_number[1])
    c_ccd_inc = dblarr(pix_number[0],pix_number[1])
    c_ccd_emi = dblarr(pix_number[0],pix_number[1])
    c_ccd_azi = dblarr(pix_number[0],pix_number[1])
    c_ccd_bri = dblarr(pix_number[0],pix_number[1])
    c_ccd_pos = dblarr(pix_number[0],pix_number[1])

    ;; 太陽直下点ベクトル
    l_ssl = dblarr(3)
    l_ssl[0] = cos(ssl_lat)*cos(ssl_lon)
    l_ssl[1] = cos(ssl_lat)*sin(ssl_lon)
    l_ssl[2] = sin(ssl_lat)

    ;; 衛星直下点
    l_ssc = dblarr(3)
    l_ssc[0] = cos(ssc_lat)*cos(ssc_lon)
    l_ssc[1] = cos(ssc_lat)*sin(ssc_lon)
    l_ssc[2] = sin(ssc_lat)

    ;; phase angle
    ;; ssl lon - ssc lon ベクトル
    tmp_l_ssl = dblarr(3)
    ;tmp_l_ssl[0] = cos(ssl_lat)*cos(ssl_lon-ssc_lon)
    ;tmp_l_ssl[1] = cos(ssl_lat)*sin(ssl_lon-ssc_lon)
    tmp_l_ssl[0] = cos(ssl_lat)*cos(ssc_lon-ssl_lon)
    tmp_l_ssl[1] = cos(ssl_lat)*sin(ssc_lon-ssl_lon)
    tmp_l_ssl[2] = sin(ssl_lat)
    ;; Rotate
    ssl_vec = r_mat##tmp_l_ssl
    ssl_vec /= sqrt(total(ssl_vec^2.))

    ;; Phase angle
    c_ccd_phase = acos(-c_ccd_Lv[*,*,0]*ssl_vec[0] $
                       -c_ccd_Lv[*,*,1]*ssl_vec[1] $
                       -c_ccd_Lv[*,*,2]*ssl_vec[2])
    ;tvscl,congrid(c_ccd_phase,200,200)
    ;stop

    R2 = Rv^2.
    D2 = Crange^2.
    DR2 = 2.*Rv*Crange
    l_ccd = dblarr(3)
    for j=0,pix_number[1]-1,1 do begin
      for i=0,pix_number[0]-1,1 do begin
        if c_ccd_Lv_ang[i,j] lt Limit_ang then begin
          c_ccd_pos[i,j] = 1.
          
          tmp_vec[0]=V_sx[i,j]
          tmp_vec[1]=V_sy[i,j]
          tmp_vec[2]=V_sz[i,j]
          latlon_vec = i_r_mat##tmp_vec
          c_ccd_lat[i,j]=asin(latlon_vec[2])
          c_ccd_lon[i,j]=-atan(latlon_vec[1],latlon_vec[0])+ssc_lon ;; "-"をつける理由未解決(2012/02/25)
          
          l_ccd[0] = cos(c_ccd_lat[i,j])*cos(c_ccd_lon[i,j])
          l_ccd[1] = cos(c_ccd_lat[i,j])*sin(c_ccd_lon[i,j])
          l_ccd[2] = sin(c_ccd_lat[i,j])

          ;; Incident Angle
          c_ccd_inc[i,j] = acos(total(l_ssl*l_ccd))

          ;; Emission Angle
          tmp_angle = acos(total(l_ssc*l_ccd))
          tmp_X = sqrt(R2+D2-DR2*cos(tmp_angle))
          tmp_angle2 = acos((R2+tmp_X^2.-D2)/(2.*Rv*tmp_X))
          c_ccd_emi[i,j] = !dpi-tmp_angle2

         
;          ;; Bright
;          imu = cos(c_ccd_inc[i,j])
;          emu = cos(c_ccd_emi[i,j])
;          if imu gt 0. then c_ccd_bri[i,j] = 0.59/!dpi*(abs(emu*imu))^0.90/emu*(1-exp(-imu/0.0547))/(1-exp(-emu/0.0039))
        endif
      endfor
    endfor
    
    ;; Azimuthal angle ;; SPモデルにはいらない
    imu = cos(c_ccd_inc)
    emu = cos(c_ccd_emi)
    inu = sin(c_ccd_inc)
    enu = sin(c_ccd_emi)
    azi_pos = where(inu*enu ne 0)
    tmp_azi = inu*0.
    tmp_azi[azi_pos] = (cos(c_ccd_phase[azi_pos])-imu[azi_pos]*emu[azi_pos])/ $
                       abs(inu[azi_pos]*enu[azi_pos])
    tmp_azi = (tmp_azi < 1 ) > (-1)
    c_ccd_azi = tmp_azi*0.
    c_ccd_azi[azi_pos] = acos(tmp_azi[azi_pos])
    
    ;; Longitudeの整形 -pi < lon < pi
    c_ccd_lon = c_ccd_lon + 2*!dpi*(c_ccd_lon lt -!dpi) - 2*!dpi*(c_ccd_lon gt !dpi)
    ;c_ccd_lon = -c_ccd_lon

    ;; 宇宙空間の値を外れ値に
    c_ccd_lon = c_ccd_lon*(c_ccd_pos ne 0) - !dpi*(c_ccd_pos eq 0)
    c_ccd_lat = c_ccd_lat*(c_ccd_pos ne 0) - !dpi*(c_ccd_pos eq 0)
    c_ccd_emi = c_ccd_emi*(c_ccd_pos ne 0) - !dpi*(c_ccd_pos eq 0)
    c_ccd_phase = c_ccd_phase*(c_ccd_pos ne 0) - !dpi*(c_ccd_pos eq 0)
    c_ccd_inc = c_ccd_inc*(c_ccd_pos ne 0) - !dpi*(c_ccd_pos eq 0)

    c_ccd_geo = dblarr(pix_number[0],pix_number[1],6)
    c_ccd_geo[*,*,0] = c_ccd_lon
    c_ccd_geo[*,*,1] = c_ccd_lat
    c_ccd_geo[*,*,2] = c_ccd_inc
    c_ccd_geo[*,*,3] = c_ccd_emi
    c_ccd_geo[*,*,4] = c_ccd_phase
    c_ccd_geo[*,*,5] = surface_dis

    error = 0.
    return,c_ccd_geo
end

;;
;; Yokota et al., 2011のアンシナリファイルに記述されている種々角度依存係数を読み込む
;;
function map_read_photometric_params, ifname,pix_n,wav,b0,db0,h,dh,c,dc,g,dg,r_mean
  tmp =''
  pix_n =0l
  tmp_wav = 0.
  tmp_B0 = 0.
  tmp_dB0 = 0.
  tmp_h = 0.
  tmp_dh = 0.
  tmp_c = 0.
  tmp_dc = 0.
  tmp_g = 0.
  tmp_dg = 0.
  tmp_mean = 0.
  
  pix_n = lonarr(181)
  wav = dblarr(181)
  b0 = dblarr(181)
  db0 = dblarr(181)
  h = dblarr(181)
  dh = dblarr(181)
  c = dblarr(181)
  dc = dblarr(181)
  g = dblarr(181)
  dg = dblarr(181)
  r_mean = dblarr(181)
  
  close,1
  openr,1,ifname
  readf,1,tmp
  count = 0l
  while(not EOF(1)) do begin
    readf,1,tmp_pix_n,tmp_wav,tmp_B0,tmp_dB0,tmp_h,tmp_dh,tmp_c,tmp_dc,tmp_g,tmp_dg,tmp_mean
    ;print,pix_n,tmp_wav,tmp_B0,tmp_h,tmp_c,tmp_mean
    
    pix_n[count] = tmp_pix_n
    wav[count] = tmp_wav
    b0[count] = tmp_b0
    db0[count] = tmp_db0
    h[count] = tmp_h
    dh[count] = tmp_dh
    c[count] = tmp_c
    dc[count] = tmp_dc
    g[count] = tmp_g
    dg[count] = tmp_dg
    r_mean[count] = tmp_mean
    count++
  endwhile
  close,1
  
  pix_n = pix_n[0:count-1]
  wav = wav[0:count-1]
  b0 = b0[0:count-1]
  db0 = db0[0:count-1]
  h = h[0:count-1]
  dh = dh[0:count-1]
  c = c[0:count-1]
  dc = dc[0:count-1]
  g = g[0:count-1]
  dg = dg[0:count-1]
  r_mean = r_mean[0:count-1]

 
  return,count
end

;;
;; Henyey-Greenstein function
;;
function map_photomertic_function,phase,B0,h,c,g
  B = B0/(1.+tan(phase/2.)/h)
  P_HGp = (1.-g^2.)/(1.+g^2.-2.*g*cos(phase))^(3./2.)
  P_HGm = (1.-g^2.)/(1.+g^2.-2.*(-g)*cos(phase))^(3./2.)
  P = (1.-c)/2.*P_HGp+(1.+c)/2.*P_HGm
  ccd_f = (1.+B)*P
  
  return, ccd_f
end

;;
;; Rcc データを読み込む
;;
function read_sp_reflectance_rcc,ifname,wav,ref_rcc
  tmp_l=''
  tmp_wav = 0.0
  tmp_ref_rcc = 0.0
  wav = dblarr(10000l)
  ref_rcc = dblarr(10000l)
  close,1
  openr,1,ifname
  readf,1,tmp_l
  count=0l
  while(not EOF(1) and count le 2500) do begin
    readf,1,tmp_wav,tmp_ref_rcc
    wav[count]=tmp_wav
    ref_rcc[count]=tmp_ref_rcc
    count++
  endwhile
  close,1
  wav=wav[0:count-1]
  ref_rcc=ref_rcc[0:count-1] ;& ref_m3 = smooth(ref_m3,20,/edge_truncate)

  return,0

end

;;
;; main
;;
pro lunar_map_plot,obs_geo,ccd_ref $
                  ,ccd_geo=ccd_geo,dtheta=dtheta,pix_number=pix_number,ccd_center=ccd_center $
                  ,datadir=datadir,correct_reflectance=correct_reflectance $
                  ;; マップを外に引継ぎ用 ;;
                  ,map_sp_ref_out=map_sp_ref_out,map_sp_ref_ori=map_sp_ref_ori

  ;; 入力がなかった時の処理:: テスト処理
  if n_elements(obs_geo) eq 0 then begin
    ssl_lat_deg = 0.
    ssl_lon_deg = 0.
    ssc_lat_deg = 0.
    ssc_lon_deg = -30.
    NA_deg = 90.
    distance =  384400 -  12756.274/2.
    phi_deg =  0/3600. ;; 天体中心と視線ベクトルのなす角
    lam_deg =  0. ;; x軸、視線中心、天体中心のなす角
  endif else begin
    ssl_lat_deg = obs_geo.ssl_latitude
    ssl_lon_deg = obs_geo.ssl_longitude
    ssc_lat_deg = obs_geo.ssc_latitude
    ssc_lon_deg = obs_geo.ssc_longitude
    NA_deg = obs_geo.N_azimuth
    distance =  obs_geo.TO_distance
    phi_deg =  obs_geo.phi ;; 天体中心と視線ベクトルのなす角
    lam_deg =  obs_geo.lam ;; x軸、視線中心、天体中心のなす角
    ;stop
  endelse
  
  if n_elements(dtheta) eq 0 then dtheta = 5./3600.*!dpi/180. ;; rad/pix CCDの出力に関係する
  if n_elements(pix_number) le 1 then pix_number = [512,512] ;; pixel数
  if n_elements(ccd_center) eq 0 then ccd_center = [255.5, 255.5]
  phot_param = dblarr(12)

  ;; 天体半径、天体までの距離
  act_Rm = 3474.3/2d ;; km
  Crange = distance

  pix2sec = 3600d*180d/!dpi*dtheta
  rad2sec = 3600d*180d/!dpi  
  sec2rad = 1d/rad2sec  

  ;;;;;;;;;;;;;;;;;;;;;;
  ;; Input file names ;;
  ;;;;;;;;;;;;;;;;;;;;;;
  if n_elements(datadir) eq 0 then begin
    fldname ='C:\work\ENVI_IDL\LunaCal\Data\'
  endif else begin
    fldname = datadir
  endelse
  
  ;; SP data 1x1
  ;ifname_ref = fldname +'avg_cube_1000s-7000s_slct_smtcont110303.img'
  ;; SP data 0.5x0.5
  ifname_ref = fldname +'avg_cube_1000s-7000s_selected_ip110225.img'
  ;; 混乱の下となると思うので元データをいじったデータは使わない。代わりにプログラムの中で更新する。
  ;ifname_ref = fldname +'avg_cube_1000s-7000s_selected_ip110225_new.img'

  ;; SP reflectance collect rcc
  ;ifname_ref_rcc = fldname +'sp_reflectance_collect.txt'

  ;; Albedo group map
  ifname_albedo = fldname + 'albedo_group_05x05.dat'
  ;; Albedo group毎の角度依存パラメータ記述ファイル
  ifname_albedo_g_h = fldname+'High_albedo_sel.txt'
  ifname_albedo_g_m = fldname+'Mid_albedo_sel.txt'
  ifname_albedo_g_l = fldname+'Low_albedo_sel.txt'

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Longitude Latitude ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;
  gpd = long(1./0.5) ;; grid per degree
  lon = dblarr(360*gpd,180*gpd)
  lat = dblarr(360*gpd,180*gpd)

  for i=0, 360*gpd-1, 1 do begin
    lon[i,*] = double(i)/gpd-180
  endfor
  for j=0, 180*gpd-1, 1 do begin
    lat[*,j] = double(j)/gpd-90
  endfor
  lon *= !dpi/180.
  lat *= !dpi/180.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Lunar Refrectance and Albedo_g ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 0.5x0.5 deg Refrectance and Albedo group (data_type 2 = 4-byte signed interger)
  map_sp_ref = read_binary(ifname_ref, data_type=5, data_dims=[360*2,180*2,160])
  map_albedo_g = read_binary(ifname_albedo, data_type=2, data_dims=[360*2,180*2])

  siz_sp_ref = size(map_sp_ref)
  map_sp_ref = shift(map_sp_ref,siz_sp_ref[1]/2.,0,0)
  map_sp_ref = reverse(map_sp_ref,2)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; SP反射率のCollection 係数を与えておく ;; redding 対策
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if keyword_set(correct_reflectance) then begin
    ;error = read_sp_reflectance_rcc(ifname_ref_rcc,wav_rcc,ref_rcc)
    stop
  endif else begin
    ;; 必要のないときはすべて1
    ref_rcc = 1.+0*dindgen(siz_sp_ref[3])
  endelse  

  ;; Albedo グループごとの波長依存位相関数パラメータの読み込み
  wav = channel_wavelength()
  count = map_read_photometric_params(ifname_albedo_g_h,pix_n,wav,h_b0,h_db0,h_h,h_dh,h_c,h_dc,h_g,h_dg,h_r_mean) 
  count = map_read_photometric_params(ifname_albedo_g_m,pix_n,wav,m_b0,m_db0,m_h,m_dh,m_c,m_dc,m_g,m_dg,m_r_mean) 
  count = map_read_photometric_params(ifname_albedo_g_l,pix_n,wav,l_b0,l_db0,l_h,l_dh,l_c,l_dc,l_g,l_dg,l_r_mean) 

  ;;;;;;;;;;;;;;;;;;;;  
  ;; Lunar Geometry ;;
  ;;;;;;;;;;;;;;;;;;;;
  ; Sub solor latitude, longitude
  ssl_lat = ssl_lat_deg*!dpi/180.
  ssl_lon = ssl_lon_deg*!dpi/180.

  ; Sub spacecraft latitude
  ssc_lat = ssc_lat_deg*!dpi/180.
  ssc_lon = ssc_lon_deg*!dpi/180.
  ; North Azimuth defined by Kouyama Phd thesis
  Na =NA_deg*!dpi/180.

  phi = phi_deg*!dpi/180.
  lam = lam_deg*!dpi/180.

  ;print,"Distance/Rvenus : ",Crange/Rv, Crange,Rv

  theta = asin(act_Rm/Crange)
  phi1 = phi-theta
  phi2 = phi+theta

;  ;; 天体中心の投影点
  xc = -tan(phi)*cos(lam)
  yc = -tan(phi)*sin(lam)

  ;; focal_lenght = 1に規格化
  fl = 1.

  ;; 楕円中心の位置 (radian)
  e_x0 = -(tan(phi1)+tan(phi2))/2.*cos(lam)
  e_y0 = -(tan(phi1)+tan(phi2))/2.*sin(lam)
  
  ;; CCD平面上での原点から楕円中心までの距離
  r = sqrt(e_x0^2.+e_y0^2.)
  ;; 楕円の傾き
  ;lam = atan(y0,x0)

  ;; 長軸の長さ
  semi_a = abs(tan(phi2)-tan(phi1))/2.
  semi_b = semi_a*sqrt(cos(phi)^2.-sin(phi)^2.*tan(theta)^2.)
  print,"### Exp radius: ", semi_a/tan(dtheta)
  print,"### Exp diameter: ", semi_a*2./tan(dtheta)
  ;print,"Eccentricity : ",sqrt(1-(semi_b/semi_a)^2.),semi_a,semi_b
  
  ;; lv vector
  l_v = dblarr(3)
  abs_lv = sqrt(xc^2.+yc^2.+1.^2.)
  l_v[0] = -xc/abs_lv
  l_v[1] = -yc/abs_lv
  l_v[2] = 1./abs_lv
  
  print,"Lv Vector : ",l_v[0],l_v[1],l_v[2] ;,sqrt(total(l_v^2.))

  ;; Venus Center
  Dis = Crange/act_Rm
  Rv = 1.
  V_c = dblarr(3)
  V_c[0] = Dis*l_v[0]
  V_c[1] = Dis*l_v[1]
  V_c[2] = Dis*l_v[2]
  
  ;; <メモ>
  ;; 天体の中心から天体の北極へ向かうベクトルを求める。
  ;; ここで使うNA(=North Azimuth)は画像上におけるReferenceベクトル(ex. x軸)と北極ベクトルのなす角ではない
  ;; 観測者基準座標系におけるベクトルの方位角。詳しいことはあかつきL3の資料やKouyama PhD. Thesis参照。
  ;; ( http://dl.dropbox.com/u/58145314/Papers/D_thesis_Kouyama_pub_main_en.pdf )
  ;; 方位角が決まると衛星直下点緯度と合わせることで仰角も決まる
  ;; (衛星直下点緯度が極域近くでは一意に決まらない時もあるが)

  ;; (2012/04/26追記) そもそもSpiceを使うなら観測者基準座標系で北極に向かうベクトルを計算して用意してきてもよい。
  ;; (2012/08/18追記) Terra/ASTERのfk, pck kernelを手に入れられていないので今回は自分で計算する
  ;;                  Terraのspk kernel はTwo line elements (tle)から自分で用意した

  alpha = atan((l_v[0]*cos(Na)+l_v[1]*sin(Na))/l_v[2])
  AA = l_v[2]
  BB = (l_v[0]*cos(Na)+l_v[1]*sin(Na))
  ;print,"Alpha, sqrt", alpha*180./pi,-sin(ssc_lat)/sqrt(AA^2.+BB^2.)
  ele = asin(-sin(ssc_lat)/sqrt(AA^2.+BB^2.))-alpha
  if abs(ele) gt !dpi/2. then begin
    ele = (!dpi-abs(ele))*ele/(abs(ele))
    Na = Na+!dpi
  endif
  tmp_al = -asin(sin(alpha)*sqrt(AA^2.+BB^2.))
  ;print,"North Vector Elevation in CCD Frame Coordinate",ele*180./pi,tmp_al*180./pi
  
  ;; Sphere basic 3 axes
  Axes = dblarr(3,3)
  ;; Z :: From Center to North Pole
  Axes[2,0] = cos(ele)*cos(Na)
  Axes[2,1] = cos(ele)*sin(Na)
  Axes[2,2] = sin(ele)
  ;; Y = (-l_v) × Z
  Axes[1,0] = 1./cos(ssc_lat)*(-l_v[1]*Axes[2,2]+l_v[2]*Axes[2,1])
  Axes[1,1] = 1./cos(ssc_lat)*(-l_v[2]*Axes[2,0]+l_v[0]*Axes[2,2])
  Axes[1,2] = 1./cos(ssc_lat)*(-l_v[0]*Axes[2,1]+l_v[1]*Axes[2,0])
  ;; X = Z × Y
  Axes[0,0] = (Axes[2,1]*Axes[1,2]-Axes[2,2]*Axes[1,1])
  Axes[0,1] = (Axes[2,2]*Axes[1,0]-Axes[2,0]*Axes[1,2])
  Axes[0,2] = (Axes[2,0]*Axes[1,1]-Axes[2,1]*Axes[1,0])

  L_N = Axes[2,*]
  
  pr = l_v[0]*Axes[0,0]+l_v[1]*Axes[0,1]+l_v[2]*Axes[0,2]
  pr = pr > (-1) & pr = pr < 1.
  
  error = latlon_map(lat,lon,Crange,act_Rm,L_v,L_N,ssl_lat,ssl_lon,ssc_lat,ssc_lon $
                    ,map_distance,map_inc,map_emi,map_phase)
  print,"###", max(map_phase)*180./!dpi,mean(map_phase)*180./!dpi,min(map_phase)*180./!dpi


  ;; 各角度のマップを出力 ;;
  ;wset,1
  ;erase
  ;xs = !D.x_size
  ;ys = !D.y_size
  ;tvscl,congrid(map_phase,ys*2./3.,ys/3.)
  ;tvscl,congrid(map_emi,ys*2./3.,ys/3.),0,ys/3.
  ;tvscl,congrid(map_inc,ys*2./3.,ys/3.),0,ys*2./3.
  ;wset,0
  
  ;;;;;;;;;; ここまで地図座標上での計算 ;;;;;;;;;;;;;;;;;

  ;;;;;;;;;; ここからCCD投影 ;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;
  ;; ccd_geo をうまく更新してやれば歪み考慮後の画像に適応できる? ;; 2014/04/04 ;;
  ;; 以下に楽に実装するか！
  ;;
    
  ;;;;;;;;;;;;;;;;;;
  ;; CCD geometry ;;
  ;;;;;;;;;;;;;;;;;;
  map_siz = size(map_sp_ref)
  
  if n_elements(ccd_geo) eq 0 then begin
    ccd_geo = map_ccd_latlon(dtheta,pix_number,ccd_center, distance,act_Rm,L_v,L_N,ssl_lat,ssl_lon,ssc_lat,ssc_lon)
  endif else begin
    ccd_siz = size(ccd_geo)
    pix_number = [ccd_siz[1],ccd_siz[2]]
  endelse

  ccd_lon = ccd_geo[*,*,0]
  ccd_lat = ccd_geo[*,*,1]
  ccd_inc = ccd_geo[*,*,2]
  
  ccd_pos = where(cos(ccd_inc) gt -1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Observed refrectance ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Interpolation用の座標用意

  ;; Yokotaさん提供マップは0.25 degから始まることによる補正 (2012/06/21)
  zero_point_correct = -0.25

  ;; "-0.25"の符号がマイナスの理由：
  ;; Yokota mapは0.25deg が 原点 (あるいは最初のpixelの中心点)
  ;; 0 deg は原点より負側に出っ張っている

  map_x = (ccd_lon[ccd_pos]*180./!dpi + zero_point_correct)*2. +180.*2
  map_x = map_x + 360.*2*(map_x lt 0)
  map_y = (ccd_lat[ccd_pos]*180./!dpi + zero_point_correct)*2. +90.*2

  ccd_ref = dblarr(pix_number[0],pix_number[1],map_siz[3])
  tmp_ccd_ref = dblarr(pix_number[0],pix_number[1])

  ;; 初期化 & 配列用意
  X_L = dblarr(map_siz[1],map_siz[2])
  map_f = dblarr(map_siz[1],map_siz[2])

  map_sp_ref_ori=dblarr(map_siz[1],map_siz[2],map_siz[3])
  map_sp_ref_out=dblarr(map_siz[1],map_siz[2],map_siz[3])

  ;; Arrays for Reference condition (30, 0, 30)
  map_f_corr = dblarr(map_siz[1],map_siz[2])
  map_phase_corr = dblarr(map_siz[1],map_siz[2])
  map_phase_corr[*] = 30.*!dpi/180. ;; 基準は30° & 配列として用意

  ;;;;;;;;;;;;;;;;;;;;;;;  
  ;; Lunar Lambert law ;;
  ;;;;;;;;;;;;;;;;;;;;;;;  
  i_deg = (map_phase*180./!dpi);*(map_inc gt -!dpi)*(map_inc le !dpi/2.)
  ;i_deg = (map_inc*180./!dpi)*(map_inc gt -!dpi)*(map_inc le !dpi/2.)
  c1 = -0.019
  c2 = 0.242d-3
  c3 = -1.46d-6
  Lf = (1.+c1*i_deg+c2*i_deg^2.+c3*i_deg^3.) > 0. ;; Limb darkening function

;  map_pos = where(((cos(map_inc)+cos(map_emi)) gt 0) and (cos(map_inc) ge 0))
;  X_L[map_pos] = 2.*Lf[map_pos]*cos(map_inc[map_pos])/(cos(map_inc[map_pos])+cos(map_emi[map_pos])) $
;               +(1.-Lf[map_pos])*cos(map_inc[map_pos])

  ;; 値域に怪しい部分があったので修正 2015.11.10 ;;
  map_pos = where(((cos(map_emi)) gt 0) and (cos(map_inc) ge 0))
  X_L[map_pos] = 2.*Lf[map_pos]/(1.+cos(map_emi[map_pos])/cos(map_inc[map_pos])) $
                 +(1.-Lf[map_pos])*cos(map_inc[map_pos])
  ;print,max(X_L),min(X_L)
  ;print,max(Lf),min(Lf)
  ;tvscl,Lf
  ;stop
  
  ;; Reference condition 30, 0, 30 ;;
  Lf_corr = (1.+c1*30.+c2*30.^2.+c3*30.^3.)
  X_L_corr = 2.*Lf_corr*cos(30.*!dpi/180.)/(cos(30.*!dpi/180.)+cos(0.))+(1.-Lf_corr)*cos(30.*!dpi/180.)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Lunar Refrectance with Phase dependency ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  pu_l = 10
  ;for i=pu_l, pu_l,1 do begin  
  ;for i=0, 80-1,1 do begin
  for i=0, map_siz[3]-1,1 do begin

    ;; アルベドグループごとの位相角依存パラメータ
    ;; 波長依存性があるので波長ごとの値を使う
    phot_param = [h_b0[i],h_h[i],h_c[i],h_g[i] $
                 ,m_b0[i],m_h[i],m_c[i],m_g[i] $
                 ,l_b0[i],l_h[i],l_c[i],l_g[i] ]

    ;; Albedo groupごとの係数をつかってfを計算
    ;; High albedo region
    h_pos = where(map_albedo_g eq 3)
    B0 = phot_param[0]
    h = phot_param[1]
    c = phot_param[2]
    g = phot_param[3]
    h_map_f = map_photomertic_function(map_phase[h_pos],B0,h,c,g)
    h_map_f_corr = map_photomertic_function(map_phase_corr[h_pos],B0,h,c,g)

    ;; Mid albedo region
    m_pos = where(map_albedo_g eq 2)
    B0 = phot_param[4]
    h = phot_param[5]
    c = phot_param[6]
    g = phot_param[7]
    m_map_f = map_photomertic_function(map_phase[m_pos],B0,h,c,g)
    m_map_f_corr = map_photomertic_function(map_phase_corr[m_pos],B0,h,c,g)

    ;; Low albedo region
    l_pos = where(map_albedo_g eq 1)
    B0 = phot_param[8]
    h = phot_param[9]
    c = phot_param[10]
    g = phot_param[11]
    l_map_f = map_photomertic_function(map_phase[l_pos],B0,h,c,g)
    l_map_f_corr = map_photomertic_function(map_phase_corr[l_pos],B0,h,c,g)

    ;; 結合
    map_f[h_pos]=h_map_f
    map_f[m_pos]=m_map_f
    map_f[l_pos]=l_map_f

    map_f_corr[h_pos]=h_map_f_corr
    map_f_corr[m_pos]=m_map_f_corr
    map_f_corr[l_pos]=l_map_f_corr
    map_f_corr[where(map_f_corr eq 0.)] = 1.

    ;; 反射率へ, 補正係数かけておく
    map_sp_ref_ori[*,*,i] = map_sp_ref[*,*,i]*ref_rcc[i]
    map_sp_ref_out[*,*,i] = (map_sp_ref[*,*,i]*ref_rcc[i])*X_L/X_L_corr*map_f/map_f_corr
    ;map_sp_ref_out[*,*,i] = (map_sp_ref[*,*,i]*ref_rcc[i]);*X_L/X_L_corr

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; CCD 画像へ Interpolation ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    tmp_map = map_sp_ref_out[*,*,i]

;    if abs(ssc_lon_deg) lt 90 then begin
;      result = interpolate(tmp_map,map_x,map_y,missing=0,/double)
;      tmp_ccd_ref[ccd_pos]=result
;
;    endif else begin
;      ;; 切れ目対策が必要 ;;
;      tmp_map_x = map_x + 360.
;      tmp_map_x = tmp_map_x - 360.*2*(tmp_map_x gt 360.*2)
;      tmp_map = shift(tmp_map,360,0)      
;      result = interpolate(tmp_map,tmp_map_x,map_y,missing=0)
;      tmp_ccd_ref[ccd_pos]=result
;    endelse

    ;; 経度0°で値がないことを防ぎたい ;;
    map_siz = size(tmp_map)
    cyc_tmp_map = dblarr(map_siz[1]+2,map_siz[2])
    ;; cyclic処理 ;;
    cyc_tmp_map[0,*]=tmp_map[map_siz[1]-1,*]
    cyc_tmp_map[1:map_siz[1],*]=tmp_map
    cyc_tmp_map[map_siz[1]+1,*]=tmp_map[0,*]
    tmp_map = cyc_tmp_map

    ;; cyclicを意識した内挿 ;;
    result = interpolate(tmp_map,map_x+1,map_y,missing=0,/double)
    tmp_ccd_ref[ccd_pos]=result

    ccd_ref[*,*,i]=tmp_ccd_ref

    if i mod 10 eq 0 then print,i,wav[i]

  endfor

  ;; Output ;;
  loadct,39,/silent
  !p.position = [0.175,0.175,0.9,0.9]

  ;; 陰影のない月面表示 sample;;
  output_moon_surface = "N"
  if output_moon_surface eq "Y" then begin
    wset,0
    tmp_map = map_sp_ref[*,*,40]

    ;; 月表側用
    ;ofname = 'C:\work\ENVI_IDL\SELENE\Products\nearside.tif'
    ;tmp_map_x = map_x

    ;; 月裏側用
    ofname = 'C:\work\ENVI_IDL\SELENE\Products\farside.tif'
    tmp_map_x = map_x + 360.
    tmp_map_x = tmp_map_x - 360.*2*(tmp_map_x gt 360.*2)
    tmp_map = shift(tmp_map,360,0)

    result = interpolate(tmp_map,tmp_map_x,map_y,missing=0)
    tmp_ccd_ref[ccd_pos]=result
    sample_ccd_ref = tmp_ccd_ref

    output = sample_ccd_ref
    output = long(output/0.2*255.)
    output = (output > 0) < 255
    loadct,0
    tvscl,reverse(output,2)
    loadct,39

    ;print,max(output),mean(output[where(output ne 0)]),wav[40]
    ;write_tiff,ofname,reverse(output,1)
    wset,0
  endif

  ;; 各種ジオメトリ用
  wset,1
  erase
  xs = !D.x_size
  ys = !D.y_size
  tvscl,congrid(map_phase,ys*2./3.,ys/3.)
  tvscl,congrid(map_emi,ys*2./3.,ys/3.),0,ys/3.
  tvscl,congrid(map_inc,ys*2./3.,ys/3.),0,ys*2./3.

  wset,2
  erase
  ;loadct,0,/silent
  xs = !D.x_size
  ys = !D.y_size
  tvscl,congrid(map_sp_ref[*,*,40]<0.2,ys*2./3.,ys/3.)
  tvscl,congrid(map_albedo_g,ys*2./3.,ys/3.),0,ys/3.
  tvscl,congrid(map_sp_ref_out[*,*,pu_l]<0.25,ys*2./3.,ys/3.),0,ys/3.*2
  wset,0

  ;;;;;;;;;;;;
  ;; 論文用 ;;
  ;;;;;;;;;;;;
  ;window,3,xs=720,ys=360
  ;loadct,0,/silent
  ;;tvscl,congrid(map_sp_ref[*,*,40]<0.3,720,360)
  ;tvscl,congrid(map_sp_ref_out[*,*,pu_l]<0.2,720,360)
  ;plots,ssc_lon_deg*2.+360,ssc_lat_deg*2.+180,/device,psym=1,symsize=3,color=255,thick=2
  ;plots,ssl_lon_deg*2.+360,ssl_lat_deg*2.+180,/device,psym=7,symsize=3,color=255,thick=2

  ;stop

;  wset,2
;  loadct,39,/silent
;  tvscl,congrid(ccd_ref[*,*,10]<0.2,300,300),0,0
;  loadct,39,/silent
;  wset,0
  return
end


