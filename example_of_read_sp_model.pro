;;
;; Evaluating Moon surface radiance based on map version of SP model
;; Coded by T. Kouyama (t.kouyama@aist.go.jp)
;;
;; If you use this code for your research purpose,
;; please cite below papers:
;; Yokota et al., 2011: Lunar photometric properties at wavelengths 0.5-1.6μm acquired by SELENE Spectral Profiler and their dependency on local albedo and latitudinal zones, Icarus, 215, 639-660
;; Ogohara et al., 2012: Automated cloud tracking system for the Akatsuki Venus Climate Orbiter data, Icarus, 217, 661-668
;; Kouyama et al., 2016: Development of an application scheme
;; for the SELENE/SP lunar reflectance model for radiometric calibration of hyperspectral and multispectral sensors, Planet. Space Sci., 124, 76-83
;;

Forward_function read_photometric_params
Forward_function read_sp_reflectance_rcc
Forward_function SP_channel_wavelength
Forward_function f_get_irradiance_integ
Forward_function SP_sol_irradiance

Forward_function photomertic_function_map
Forward_function geometry_map
Forward_function radiance_factor_map


;;
;; main
;;
pro example_of_read_SP_model

  ;; set plot window
  window,0,xs=360*2,ys=180*3

  ;;;;;;;;;;;;;;;;;;;;;;
  ;; Input file names ;;
  ;;;;;;;;;;;;;;;;;;;;;;
  ;; Please replace below path with your directory name where SP files are contained.
  ;ifldname ='C:\work\ENVI_IDL\sample\SP_model_example\'
  ifldname = './'
  ofldname = './'

  ;; SP map data with 0.5x0.5 grid interval (double precision, idl data_type=5) ;;
  ifname_ref = ifldname +'avg_cube_1000s-7000s_selected_ip110225.img'

  ;; Note: Definition of SP map coordinate  (Lon, Lat)
  ;; (0.25, 89.75), (0.75, 89.75), ... (359.75,89.75) ;; North pole
  ;; (0.25, 89.25), (0.75, 89.25), ... (359.75,89.25)
  ;; ...
  ;; (0.25, 0.25), (0.75, 0.25), ... (359.75,0.25)
  ;; (0.25, -0.25), (0.75, -0.25), ... (359.75,-0.25)
  ;; ...
  ;; (0.25, -89.25), (0.75, -89.25), ... (359.75,-89.25)
  ;; (0.25, -89.75), (0.75, -89.75), ... (359.75,-89.75) ;; South pole

  ;; Albedo group map (this file is made by Kouyama based on Yokota's definition in Yokota et al 2011) ;;
  ifname_albedo = ifldname + 'albedo_group_05x05.dat'

  ;; Parameters of i, e, alpha dependence for each albedo group ;;
  ifname_albedo_g_h = ifldname+'High_albedo_sel.txt'
  ifname_albedo_g_m = ifldname+'Mid_albedo_sel.txt'
  ifname_albedo_g_l = ifldname+'Low_albedo_sel.txt'

  ;; Solar irradiance data from Gueymard model (used for generating SP model)
  solfname = ifldname+'Gueymard.txt'

  ;;
  ;; Obsevation geometry Parameters
  ;; Please change below parameters according to observation geometry of your satellite.
  ;; The observation geometry can be obtained from SPICE toolkit, for example.
  ;;

  ;; Sub-solar latitude and longitude (deg)
  ssl_lat_deg = 3.
  ssl_lon_deg = 10.

  ;; Sub-spacecraft latitude and longitude (deg)
  ssc_lat_deg = 3.
  ssc_lon_deg = 0.

  ;; Distance between observer and the Moon center (km)
  LO_distance =  384400d

  ;; Distance between the Sun and the Moon (AU)
  SL_distance_au = 1d

  print, "Assumed observation geometry: "
  print, "  Sub solar latitude: ", ssl_lat_deg
  print, "  Sub solar longitude: ", ssl_lon_deg
  print, "  Sub spacecraft latitude: ", ssc_lat_deg
  print, "  Sub spacecraft longitude: ", ssc_lon_deg
  print, "  Sun-Moon Distance: ", SL_distance_au
  print, "  Moon-Observer Distance: ",LO_distance, " [km]"


  ;;
  ;; Other parameters (basically same for any satellite)
  ;;

  ;; Number of channels of SP model ;;
  channel_n = 160l
  ;; SP's wavlength ;;
  wav = SP_channel_wavelength()
  ;; Solar irradiance
  sol_irrad = SP_sol_irradiance(solfname, wav) ;; W/m2/nm
  sol_irrad = sol_irrad*1000. ;; W/m2/um

  ;; Moon radius
  Rmoon = 3474.3/2d ;; km

  ;; North direction of the Moon in a image
  ;; 0 deg: North = X direction of a image
  ;; 90 deg: North = Y direction
  NA_deg = 90.
  Na =NA_deg*!dpi/180.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Lunar Reflectance and Albedo_g ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 0.5x0.5 deg Reflectance map ;; Double precision float
  map_sp_ref = read_binary(ifname_ref, data_type=5, data_dims=[360*2,180*2,channel_n])  
  siz_sp_ref = size(map_sp_ref)

  ;; Shift with -180 deg for visualizing ;;
  map_sp_ref = shift(map_sp_ref,siz_sp_ref[1]/2.,0,0)
  ;; North is up for IDL visualization;;
  map_sp_ref = reverse(map_sp_ref,2)

  ;; Albedo group (north is up), 2 byte integer
  map_albedo_g = read_binary(ifname_albedo, data_type=2, data_dims=[360*2,180*2])

  ;; Longitude Latitude grid ;;
  gpd = long(1./0.5) ;; grid per degree
  lon = dblarr(360l*gpd,180l*gpd)
  lat = dblarr(360l*gpd,180l*gpd)

  for i=0, 360l*gpd-1, 1 do begin
    lon[i,*] = (double(i)/gpd + 0.25)-180d
  endfor

  for j=0, 180l*gpd-1, 1 do begin
    lat[*,j] = double(j)/gpd - 90d + 0.25
  endfor
  lon *= !dpi/180.
  lat *= !dpi/180.

  ;;
  ;; Map check ;;
  ;;
  erase
  ;; Reflectance map
  tvscl,map_sp_ref[*,*,0] < 0.1, 0, 90, /device
  stop
  ;; Albedo group map
  tvscl,map_albedo_g, 0, 90, /device
  stop


  ;;  
  ;; Photometric parameters
  ;;
  ;; Paramters for High albedo ;;
  h_count = read_photometric_params(ifname_albedo_g_h,pix_n,wav,h_b0,h_db0,h_h,h_dh,h_c,h_dc,h_g,h_dg,h_r_mean)

  ;; Paramters for middle albedo ;;
  m_count = read_photometric_params(ifname_albedo_g_m,pix_n,wav,m_b0,m_db0,m_h,m_dh,m_c,m_dc,m_g,m_dg,m_r_mean)

  ;; Paramters for low albedo ;;
  l_count = read_photometric_params(ifname_albedo_g_l,pix_n,wav,l_b0,l_db0,l_h,l_dh,l_c,l_dc,l_g,l_dg,l_r_mean)

  ;; Into a one array ;;
  phot_param = dblarr(12, channel_n)
  for i=0, channel_n-1, 1 do begin
    phot_param[*,i] = [h_b0[i],h_h[i],h_c[i],h_g[i] $
                      ,m_b0[i],m_h[i],m_c[i],m_g[i] $
                      ,l_b0[i],l_h[i],l_c[i],l_g[i] ]
  endfor

  ;;;;;;;;;;;;;;;;;;;;
  ;; Lunar Geometry ;;
  ;;;;;;;;;;;;;;;;;;;;
  ; Sub solor latitude, longitude
  ssl_lat = ssl_lat_deg*!dpi/180.
  ssl_lon = ssl_lon_deg*!dpi/180.

  ; Sub spacecraft latitude
  ssc_lat = ssc_lat_deg*!dpi/180.
  ssc_lon = ssc_lon_deg*!dpi/180.

  ;; Moon center position in a image frame
  xc = 0d
  yc = 0d
  ;; (Assuming Moon center is located at center of field of view)

  ;; Or you can choose any observation condition in which Moon center is located arbital position expressed with
  ;; phi: angle between boresight vector and vector from a satellite to Moon center
  ;: lam: Rotation angle of a satellite
  ;; The Moon shape may distort a little bit according to the Moon center location in a detector
  ;; when the satellite position is nearby the Moon.

  ;phi_deg =  0d
  ;lam_deg =  0d
  ;phi = phi_deg*!dpi/180.
  ;lam = lam_deg*!dpi/180.
  ;xc = -tan(phi)*cos(lam)
  ;yc = -tan(phi)*sin(lam)
  
  ;;
  ;; Vector from a satellite to Moon center (L_v) in camera coordinate ;;
  ;;
  l_v = dblarr(3)
  abs_lv = sqrt(xc^2.+yc^2.+1.^2.)
  l_v[0] = -xc/abs_lv
  l_v[1] = -yc/abs_lv
  l_v[2] = 1./abs_lv

  ;;
  ;; North vector (L_N) of the Moon in camera coordinate.
  ;; See Ogohara et al., 2012
  ;;

  ;; North direction of the Moon in a image frame should be considered
  ;; when Moon center location is not at the image center.
  alpha = atan((l_v[0]*cos(Na)+l_v[1]*sin(Na))/l_v[2])
  AA = l_v[2]
  BB = (l_v[0]*cos(Na)+l_v[1]*sin(Na))

  ele = asin(-sin(ssc_lat)/sqrt(AA^2.+BB^2.))-alpha

  if abs(ele) gt !dpi/2. then begin
    ele = (!dpi-abs(ele))*ele/(abs(ele))
    Na = Na+!dpi
  endif

  L_N = dblarr(3)
  L_N[0] = cos(ele)*cos(Na)
  L_N[1] = cos(ele)*sin(Na)
  L_N[2] = sin(ele)

  print," ## Lv Vector : ",l_v[0],l_v[1],l_v[2]
  print," ## Ln Vector : ",l_N[0],l_N[1],l_N[2]

  ;;
  ;; By using L_v and L_n vectors,
  ;; Incident, emission, and phase angles for each grid can be estimated
  ;; These angles may be useful for specifying a location of Moon surface which is suit for some observation conditions ;;
  ;;
  print,"Inc, emi, and phase angles..."
  error = geometry_map(lat,lon,LO_distance,Rmoon,L_v,L_N $
                      ,ssl_lat,ssl_lon,ssc_lat,ssc_lon $
                      ,map_distance,map_inc,map_emi,map_phase)

  ;;
  ;; QL example ;;
  ;;
  window,1,xs=720,ys=360*3
  tvscl,map_inc * (map_emi le 90*!dpi/180.),0,360*2
  tvscl,map_emi * (map_emi le 90*!dpi/180.),0,360*1
  tvscl,map_phase * (map_emi le 90*!dpi/180.) > min(map_phase),0,0
  wset,0
  ;stop

  ;; Radiance factor at each grid ;;
  print,"Radiance factor..."
  map_radiance_factor = radiance_factor_map(map_sp_ref, map_albedo_g, phot_param $
                                             ,map_inc ,map_emi ,map_phase)

  ;; Converting radiance ;;
  map_Moon_radiance = dblarr(360*2l,180*2l,channel_n)

  ;; Quick look interval
  ql_interval = 10l

  window,2,xs=720,ys=360
  for i=0, channel_n-1, 1 do begin
    map_Moon_radiance[*,*,i] = map_radiance_factor[*,*,i] $
                             * sol_irrad[i]/SL_distance_au^2./!dpi ; W/m2/um/str

    ;; for ql ;;
    if n_elements(ql_interval) ne 0 then begin
      if i mod ql_interval eq 0 then begin
        tmp = map_Moon_radiance[*,*,i]
        tmp_pos = where(tmp gt 0)
        tmp_mean = mean(tmp[tmp_pos])

        erase
        tvscl, tmp < tmp_mean * 2d,/device
        wait,0.1
        print,"Wavelength, Mean radiance: ", wav[i], tmp_mean
        ;stop

      endif
    endif

  endfor
 
  wset,0
  ;; Plot an example of radiance at (lon, lat) = (-0.25, -0.25) ;;
  plot,wav, map_Moon_radiance[359,179,*] $
      ,xtitle="Wavelength (nm)", ytitle="Radiance (W/m2/sr/um)" $
      ,title = "Radiance at Lat = -0.25, Lon = -0.25"

  ;;
  ;; Output exapmle: Mutli-band Tiff (Float)
  ;; Float tiff can be opened with, for example, QGIS.
  ;;
  ofname_result_tiff = ofldname + 'radiance_map_SP.tiff'
  write_tiff,ofname_result_tiff, transpose(map_Moon_radiance,[2,0,1]), /float

  return
end



;;
;; Read photometric parameters in supporting information of Yokota et al., 2011
;;
function read_photometric_params, ifname,pix_n,wav,b0,db0,h,dh,c,dc,g,dg,r_mean
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

  n_lines = 160l
  pix_n = lonarr(n_lines)
  wav = dblarr(n_lines)
  b0 = dblarr(n_lines)
  db0 = dblarr(n_lines)
  h = dblarr(n_lines)
  dh = dblarr(n_lines)
  c = dblarr(n_lines)
  dc = dblarr(n_lines)
  g = dblarr(n_lines)
  dg = dblarr(n_lines)
  r_mean = dblarr(n_lines)

  close,1
  openr,1,ifname
  readf,1,tmp
  count = 0l
  while(not EOF(1)) do begin
    readf,1,tmp_pix_n,tmp_wav,tmp_B0,tmp_dB0,tmp_h,tmp_dh,tmp_c,tmp_dc,tmp_g,tmp_dg,tmp_mean

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

  return,count
end

;;
;; Read Rcc data
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
;; SP model data
;;
function SP_channel_wavelength,arr_size=arr_size

  wavelength =  [ $
    512.599976,  518.400024,  524.700012,  530.400024,  536.500000,  542.799988, $
    548.700012,  554.500000,  560.500000,  566.700012,  572.599976,  578.500000, $
    584.500000,  590.599976,  596.700012,  602.500000,  608.599976,  614.599976, $
    620.500000,  626.700012,  632.700012,  638.599976,  644.599976,  650.599976, $
    656.599976,  662.599976,  668.799988,  674.700012,  680.599976,  686.700012, $
    692.599976,  698.599976,  704.700012,  710.799988,  716.700012,  722.700012, $
    728.700012,  734.700012,  740.700012,  746.799988,  752.799988,  758.700012, $
    764.799988,  770.700012,  776.700012,  782.700012,  788.799988,  794.700012, $
    800.700012,  806.799988,  812.700012,  818.700012,  824.799988,  830.799988, $
    836.799988,  842.799988,  848.799988,  854.599976,  860.700012,  866.700012, $
    872.700012,  878.700012,  884.599976,  890.700012,  896.599976,  902.700012, $
    908.700012,  914.599976,  920.599976,  926.599976,  932.599976,  938.599976, $
    944.599976,  955.400024,  963.500000,  971.400024,  979.700012,  987.599976, $
    993.700012, 1003.599976, 1013.099976, 1019.500000, 1027.699951, 1035.500000, $
    1043.599976, 1051.699951, 1059.699951, 1067.800049, 1075.800049, 1083.599976, $
    1091.800049, 1099.699951, 1107.699951, 1115.900024, 1123.800049, 1131.800049, $
    1139.699951, 1147.800049, 1155.699951, 1163.800049, 1171.800049, 1179.800049, $
    1187.800049, 1195.800049, 1203.900024, 1211.900024, 1219.800049, 1227.900024, $
    1235.900024, 1244.000000, 1252.000000, 1259.800049, 1267.800049, 1275.900024, $
    1284.199951, 1292.000000, 1299.800049, 1307.800049, 1315.900024, 1323.800049, $
    1331.800049, 1339.800049, 1347.800049, 1355.800049, 1363.800049, 1371.800049, $
    1379.800049, 1387.800049, 1395.900024, 1403.800049, 1411.800049, 1419.800049, $
    1427.900024, 1435.699951, 1443.800049, 1451.900024, 1459.800049, 1467.800049, $
    1475.800049, 1483.900024, 1491.800049, 1499.800049, 1507.800049, 1515.699951, $
    1523.800049, 1531.699951, 1539.699951, 1547.699951, 1555.500000, 1563.699951, $
    1571.699951, 1579.599976, 1587.699951, 1595.699951, 1603.699951, 1611.699951, $
    1620.099976, 1628.099976, 1636.099976, 1644.199951 $
    ]

  c_size=size(wavelength)
  arr_size=c_size[1]
  return, wavelength
end

;;
;; Converting irradiance interval from input to output using trapezoidal integration
;;
function f_get_irradiance_integ,input_wav,input_rad, output_wav

  input_wav_p = shift(input_wav,-1)
  input_wav_m = shift(input_wav,1)

  obs_int_p = (input_wav_p-input_wav)
  obs_int_m = (input_wav-input_wav_m)
  obs_int = obs_int_p + obs_int_m
  n_data = n_elements(input_wav)

  siz = size(output_wav)

  ;; Work array ;;
  tmp_output_wav = dblarr(siz[1]+1)
  tmp_output_wav[0:siz[1]-1]=output_wav
  tmp_output_wav[siz[1]] = (output_wav[siz[1]-1]-output_wav[siz[1]-2])+output_wav[siz[1]-1]

  sel_rad = dblarr(siz[1])
  position = 0l
  position1 = 0l

  ;; wav 0------wav 1------wav 2
  ;;        |<--range-->|
  for i=0, siz[1]-1,1 do begin
    if i eq 0 then begin
      while(input_wav[position] le tmp_output_wav[i]) do begin
        if position ge n_data-1 then break
        position++
      endwhile

      position1 = position
      while(input_wav[position1] le (tmp_output_wav[i]+tmp_output_wav[i+1])/2.) do begin
        if position1 ge n_data-1 then break
        position1++
      endwhile
    endif else begin
      while(input_wav[position] le (tmp_output_wav[i]+tmp_output_wav[i-1])/2.) do begin
        if position ge n_data-1 then break
        position++
      endwhile

      position1 = position
      while(input_wav[position1] le (tmp_output_wav[i]+tmp_output_wav[i+1])/2.) do begin
        if position1 ge n_data-1 then break
        position1++
      endwhile
    endelse

    if position ge n_data then return, sel_rad
    ;; Simple mean ;;
    ;sel_rad[i]=mean(input_rad[position:position1])

    ;; Averaging with trapezoidal integration ;;
    if position1 eq position then begin
      if position1 eq 0 then begin
        sel_rad[i]=0d
      endif else if position eq n_data-1 then begin
        sel_rad[i]=0d
      endif else begin
        sel_rad[i]=mean(input_rad[position:position1])
      endelse
    endif else if position1 gt position then begin
      sel_rad[i]=0d
      tmp_total_wav = 0d

      for j = position, position1-1, 1 do begin
        d_wav = (input_wav[j+1]-input_wav[j])/1000d
        sel_rad[i]+=(input_rad[j+1]+input_rad[j])*d_wav/2d
        tmp_total_wav += d_wav
      endfor
      if tmp_total_wav eq 0 then begin
        print, i,position,position1
        stop
      endif
      sel_rad[i]/= tmp_total_wav

    endif

  endfor

  return,sel_rad
end

;;
;; Read solar irradiance data
;;
function SP_sol_irradiance, solfname, wav

  count=0l
  tmp_l=''
  tmp_wav = 0.0
  tmp_rad = 0.0
  n_data = 10000l
  obs_rad = fltarr(n_data)
  obs_wav = fltarr(n_data)

  openr,unit,solfname,/get_lun
  readf,unit,tmp_l
  while(not EOF(unit)) do begin
    readf,unit,tmp_wav,tmp_rad
    obs_wav[count] = tmp_wav
    obs_rad[count] = tmp_rad
    count++
  endwhile
  close,unit
  free_lun,unit

  obs_wav = obs_wav[0:count-1]
  obs_rad = obs_rad[0:count-1]

  ;; 台形積分して太陽光スペクトルを得る ;;
  sol_rad=f_get_irradiance_integ(obs_wav,obs_rad,wav) ; W/m2/nm

  return, sol_rad
end

;;
;; Henyey-Greenstein function
;;
function photomertic_function_map,phase,B0,h,c,g
  B = B0/(1.+tan(phase/2.)/h)
  P_HGp = (1.-g^2.)/(1.+g^2.-2.*g*cos(phase))^(3./2.)
  P_HGm = (1.-g^2.)/(1.+g^2.-2.*(-g)*cos(phase))^(3./2.)
  P = (1.-c)/2.*P_HGp+(1.+c)/2.*P_HGm
  ccd_f = (1.+B)*P

  return, ccd_f
end

;;
;; Geometry condition at each grid
;;
function geometry_map,lat,lon,LO_distance,Rv $
                       ,f_Lv,f_Nv $
                       ,ssl_lat,ssl_lon $
                       ,ssc_lat,ssc_lon $
                       ;; Output ;;
                       ,c_map_distance,c_map_inc,c_map_emi,c_map_phase

  ;; Sub solar position ;;
  l_ssl = dblarr(3)
  l_ssl[0] = cos(ssl_lat)*cos(ssl_lon)
  l_ssl[1] = cos(ssl_lat)*sin(ssl_lon)
  l_ssl[2] = sin(ssl_lat)

  ;; Sub spacecraft position ;;
  l_ssc = dblarr(3)
  l_ssc[0] = cos(ssc_lat)*cos(ssc_lon)
  l_ssc[1] = cos(ssc_lat)*sin(ssc_lon)
  l_ssc[2] = sin(ssc_lat)

  ;; Normalized position at each grid ;;
  l_map_x = cos(lat)*cos(lon)
  l_map_y = cos(lat)*sin(lon)
  l_map_z = sin(lat)

  ;;
  ;; Incident Angle at each grid ;;
  ;;
  tmp_prod_inc = l_ssl[0]*l_map_x+l_ssl[1]*l_map_y+l_ssl[2]*l_map_z

  c_map_inc = acos(tmp_prod_inc *(abs(tmp_prod_inc) le 1.)) * (abs(tmp_prod_inc) le 1.)

  ;;
  ;; Emission Angle at each grid ;;
  ;;
  tmp_prod_emi = l_ssc[0]*l_map_x+l_ssc[1]*l_map_y+l_ssc[2]*l_map_z
  tmp_angle = acos((tmp_prod_emi < 1 ) > (-1))

  R2 = Rv^2.
  D2 = LO_distance^2.
  DR2 = 2.*Rv*LO_distance

  tmp_X = sqrt((R2+D2-DR2*cos(tmp_angle)) >0.)
  tmp_angle2 = acos( (((R2+tmp_X^2.-D2)/(2.*Rv*tmp_X)) < 1) > (-1) )

  c_map_emi = !dpi-tmp_angle2

  ;;
  ;; Phase angle at each grid
  ;;

  ;; Moon center direction ;;
  V_center = f_Lv*LO_distance

  ;; Moon's three basis vectors seen from a satellite
  ;; please see https://www.dropbox.com/s/4duw8p6qohq36yi/Section_2-1_Kouyama_PhDThesis_rev.pdf?dl=0
  f_E1 = crossp(f_Lv, f_Nv)/cos(ssc_lat)
  f_E0 = crossp(f_E1, f_Nv)
  f_E2 = f_Nv

  ;;
  ;; Vector from a satellite to Moon surface position of each grid ;;
  ;;
  map_Lv_x = Rv*(l_map_x*f_E0[0] + l_map_y*f_E1[0] + l_map_z*f_E2[0]) + V_center[0]
  map_Lv_y = Rv*(l_map_x*f_E0[1] + l_map_y*f_E1[1] + l_map_z*f_E2[1]) + V_center[1]
  map_Lv_z = Rv*(l_map_x*f_E0[2] + l_map_y*f_E1[2] + l_map_z*f_E2[2]) + V_center[2]

  ;; Normalize ;;
  c_map_distance = sqrt(map_Lv_x^2.+map_Lv_y^2.+map_Lv_z^2.)
  map_Lv_x /= c_map_distance
  map_Lv_y /= c_map_distance
  map_Lv_z /= c_map_distance

  ;; Sun direction, sub-spacecraft longitude as 0 deg
  tmp_l_ssl = dblarr(3)
  tmp_l_ssl[0] = cos(ssl_lat)*cos(ssl_lon-ssc_lon)
  tmp_l_ssl[1] = cos(ssl_lat)*sin(ssl_lon-ssc_lon)
  tmp_l_ssl[2] = sin(ssl_lat)

  ;; Rotate
  r_mat = dblarr(3,3)
  r_mat[0,*]=f_E0
  r_mat[1,*]=f_E1
  r_mat[2,*]=f_E2

  ssl_vec = r_mat##tmp_l_ssl
  ;ssl_vec /= sqrt(total(ssl_vec^2.))

  ;; Phase angle at each grid ;;
  c_map_phase = acos(-map_Lv_x*ssl_vec[0] $
    -map_Lv_y*ssl_vec[1] $
    -map_Lv_z*ssl_vec[2])

  return,0
end

;;
;;  Converting reflectance value to radiance factor
;;
function radiance_factor_map, map_sp_ref, map_albedo_g, phot_param $
                               ;; Outputs
                               , map_inc, map_emi, map_phase

  ;; 初期化 & 配列用意
  map_siz = size(map_sp_ref)
  X_L = dblarr(map_siz[1],map_siz[2])
  map_f = dblarr(map_siz[1],map_siz[2])

  map_sp_ref_ori=dblarr(map_siz[1],map_siz[2],map_siz[3])
  map_sp_ref_out=dblarr(map_siz[1],map_siz[2],map_siz[3])

  ;; Arrays for Reference condition (30, 0, 30)
  map_f_corr = dblarr(map_siz[1],map_siz[2])
  map_phase_corr = dblarr(map_siz[1],map_siz[2])
  ;; reference alpha is  30 deg
  map_phase_corr[*] = 30.*!dpi/180.

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Lunar Lambert law ;;
  ;;;;;;;;;;;;;;;;;;;;;;;
  i_deg = (map_phase*180./!dpi)

  ;; Limb darkening function
  c1 = -0.019
  c2 = 0.242d-3
  c3 = -1.46d-6
  Lf = (1.+c1*i_deg+c2*i_deg^2.+c3*i_deg^3.) > 0.

  map_pos = where(((cos(map_emi)) gt 0) and (cos(map_inc) ge 0))

  X_L[map_pos] = 2.*Lf[map_pos]/(1.+cos(map_emi[map_pos])/cos(map_inc[map_pos])) $
    +(1.-Lf[map_pos])*cos(map_inc[map_pos])

  ;; Reference condition: (i, e, alpha) = (30, 0, 30) ;;
  Lf_corr = (1.+c1*30.+c2*30.^2.+c3*30.^3.)
  X_L_corr = 2.*Lf_corr*cos(30.*!dpi/180.)/(cos(30.*!dpi/180.)+cos(0.))+(1.-Lf_corr)*cos(30.*!dpi/180.)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Lunar Refrectance with Phase dependency ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;ql_interval = 10l

  for i=0, map_siz[3]-1,1 do begin
    tmp_phot_param = phot_param[*,i]

    ;; Albedo groupごとの係数をつかってfを計算
    ;; High albedo region
    h_pos = where(map_albedo_g eq 3)
    B0 = tmp_phot_param[0]
    h = tmp_phot_param[1]
    c = tmp_phot_param[2]
    g = tmp_phot_param[3]
    h_map_f = photomertic_function_map(map_phase[h_pos],B0,h,c,g)
    h_map_f_corr = photomertic_function_map(map_phase_corr[h_pos],B0,h,c,g)

    ;; Mid albedo region
    m_pos = where(map_albedo_g eq 2)
    B0 = tmp_phot_param[4]
    h = tmp_phot_param[5]
    c = tmp_phot_param[6]
    g = tmp_phot_param[7]
    m_map_f = photomertic_function_map(map_phase[m_pos],B0,h,c,g)
    m_map_f_corr = photomertic_function_map(map_phase_corr[m_pos],B0,h,c,g)

    ;; Low albedo region
    l_pos = where(map_albedo_g eq 1)
    B0 = tmp_phot_param[8]
    h = tmp_phot_param[9]
    c = tmp_phot_param[10]
    g = tmp_phot_param[11]
    l_map_f = photomertic_function_map(map_phase[l_pos],B0,h,c,g)
    l_map_f_corr = photomertic_function_map(map_phase_corr[l_pos],B0,h,c,g)

    ;; 結合
    map_f[h_pos]=h_map_f
    map_f[m_pos]=m_map_f
    map_f[l_pos]=l_map_f

    map_f_corr[h_pos]=h_map_f_corr
    map_f_corr[m_pos]=m_map_f_corr
    map_f_corr[l_pos]=l_map_f_corr
    map_f_corr[where(map_f_corr eq 0.)] = 1.
    

    ;; Reflectance to radiance factor
    map_sp_ref_out[*,*,i] = map_sp_ref[*,*,i]*X_L/X_L_corr*map_f/map_f_corr

    ;; for output ;;
    if n_elements(ql_interval) ne 0 then begin

      if i mod ql_interval eq 0 then begin
        print,i
        tvscl,map_sp_ref_out[*,*,i] < 0.5
        wait,0.01
      endif

    endif

  endfor


  return,map_sp_ref_out
end
