;;
;; test program for SP model simulation
;; Toru Kouyama (AIST)
;; t.kouyama@aist.go.jp
;;

;;
;; Required IDL codes should be located in the same folder
;;
@channel_wavelength.pro
@get_radiance_lism.pro
@band_response.pro
@lunar_map_plot.pro
@channel_contrib.pro
@f_get_radiance_integ.pro
@band_response.pro
@read_sp_model_bilinear_for_pub.pro

;;
;; main
;;
pro main_SP_model
  ;;
  ;; Directory names for input parameter files and outputs
  ;;
  datadir = './parameters/'
  outdir = './outputs/'

  ;;
  ;; provide observation geometry
  ;;

  AU = 149597870700d0/1000d0 ;; 1 AU (km)
  mean_MD = 384400d ;; Mean Moon-Earth distance (km)

  ;; Sun - Moon distance
  SL_distance = AU
  
  ;; Moon - Observer distance
  LO_distance = mean_MD

  ;; simulation condition example 1
  ;sub_solar_lon_deg = -5.98037d
  ;sub_solar_lat_deg = 0d
  ;sub_sc_lon_deg = -12.0d
  ;sub_sc_lat_deg = -8.0d

  ;; simulation condition example 2
  sub_solar_lon_deg = 24d
  sub_solar_lat_deg = 0d
  sub_sc_lon_deg = -6.0d
  sub_sc_lat_deg = 0d

  ;; Phase angle ;;
  tmp_obs = [cos(sub_sc_lat_deg*!dpi/180.)*cos(sub_sc_lon_deg*!dpi/180.) $
            ,cos(sub_sc_lat_deg*!dpi/180.)*sin(sub_sc_lon_deg*!dpi/180.) $
            ,sin(sub_sc_lat_deg*!dpi/180.)]

  tmp_sol = [cos(sub_solar_lat_deg*!dpi/180.)*cos(sub_solar_lon_deg*!dpi/180.) $
            ,cos(sub_solar_lat_deg*!dpi/180.)*sin(sub_solar_lon_deg*!dpi/180.) $
            ,sin(sub_solar_lat_deg*!dpi/180.)]
  phase_angle_deg = acos(total(tmp_obs*tmp_sol))*180./!dpi
  print,"Phase angle: ",phase_angle_deg
  ;stop

  ;;
  ;; NA_deg means North-pole Azimuthal angle,
  ;; meaning the angle between X direction and North pole direction from Moon center in the image frame.
  ;; -90 deg results North Pole is up.
  ;; 0 deg results North Pole direction is X direction
  ;;
  NA_deg = -90d ;; General case

  ;; Belows are always 0 for lunar calibration case ;;
  phi_deg = 0d
  lam_deg = 0d

  ;;
  ;; Create a structure
  ;;
  tags = ['ST_distance','TO_distance','ssl_longitude','ssl_latitude' $
    ,'ssc_longitude', 'ssc_latitude', 'N_azimuth', 'phi', 'lam','phase']
  obs_geo = create_struct(tags $
    ,SL_distance $
    ,LO_distance $
    ,sub_solar_lon_deg $
    ,sub_solar_lat_deg $
    ,sub_sc_lon_deg $
    ,sub_sc_lat_deg $
    ,NA_deg $
    ,phi_deg $
    ,lam_deg $
    ,phase_angle_deg $
    )


  ;;
  ;; Run simulation
  ;;
  read_sp_model_bilinear_for_pub,obs_geo $
        , out_wav, out_hyper_image, out_irad, datadir = datadir

  ;;
  ;; output disk integrated irradiance ;;
  ;;
  ;ofname_csv = outdir+'Lunar_irradiance.csv'
  ;write_csv,ofname_csv,out_wav,out_irad ,header = ['Wavelength','Lunar irradiance']

  ;; for GDL interface
  ofname_csv = outdir+'Lunar_irradiance.csv'
  openw,1,ofname_csv
  printf,1,'Wavelength, Lunar_irradiance (W/m2/um)'

  for i=0, 160-1, 1 do begin
     tmp_l = string(double(out_wav[i]))+','+string(double(out_irad[i]))
     print,tmp_l
     printf,1,tmp_l
  endfor
  close,1

  ;;
  ;; output disk resolved simulation image (hyper) with tiff or binary format
  ;; out_hyper_image = (X, Y, Channel)
  ;; note:: output file size ~ 800 MB / cube
  ;;

  ofname_base = 'simulation_image_hyper'
  
  ;;
  ;; for IDL interface
  ;;

  ;ofname_tiff = outdir+ofname_base+'.tiff'
  ;write_tiff,ofname_tiff,out_hyper_image,/double

  ;; this tiff can be read as ::
  ;;   Moon_sim_image = read_tiff(file_name)
  ;;   help, Moon_sim_image
  ;;   >  DOUBLE    = Array[800, 800, 160]
  
  ;;
  ;; for GDL interface
  ;;
  ofname_bin = outdir+ofname_base+'.bin'
  openw,1,ofname_bin
  writeu,1,out_hyper_image
  close,1

  ;; this binary file can be opened with, for example, ImageJ via reading "Raw file"
  ;; 64-bit real, 800 x 800 size, 160 images, and Little endian

  ;;
  ;; output observation settings
  ;;
  ofname_struct = outdir+ ofname_base + '_setting.txt'
  help,obs_geo,/structure,output=str_obs_geo
  openw,1,ofname_struct

  out_size = size(out_hyper_image)
  printf,1,"Output binary file information: "
  printf,1,'64 bits/pix, Little Endian '
  printf,1,"Data size (X, Y, Channel)"
  printf,1,"X  "+string(out_size[1])
  printf,1,"Y  "+string(out_size[2])
  printf,1,"Channel  "+string(out_size[3])

  printf,1,""
  for i=0, n_elements(str_obs_geo)-1, 1 do begin
    printf,1,str_obs_geo[i]    
  endfor

  
  close,1

  ;;
  ;; QL file
  ;; 
  ofname_jpeg = outdir+ ofname_base + '_ql.jpg'
  tmp_image = out_hyper_image[*,*,40]
  byte_tmp_image = byte( tmp_image / (max(tmp_image) < 110) * 255d)
  write_jpeg,ofname_jpeg, byte_tmp_image

  return
end
