;;
;; test program for SP model simulation test
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
  ;; Directory names for input parameter files and output
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
  ;sub_solar_lon_deg = -5.98037
  ;sub_solar_lat_deg = 0d
  ;sub_sc_lon_deg = -12.0
  ;sub_sc_lat_deg = -8.0

  ;; simulation condition example 2
  sub_solar_lon_deg = 24d
  sub_solar_lat_deg = 0d
  sub_sc_lon_deg = -6.0
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

  ;; direction from Moon center to North pole in the image frame ;;
  NA_deg = -90;; General case

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
  ofname_csv = outdir+'sample_correct.csv'
  write_csv,ofname_csv,out_wav,out_irad ,header = ['Wavelength','Lunar irradiance']

  ;;
  ;; output disk resolved simulation image (hyper) with tiff format
  ;; note:: output file size ~ 800 MB / cube
  ;;
  ofname_tiff = outdir+'simulation_image_hyper.tiff'
  ;; out_hyper_image = (X, Y, Channel)
  write_tiff,ofname_tiff,out_hyper_image,/double

  ;; this tiff can be read as ::
  ;;   Moon_sim_image = read_tiff(file_name)
  ;;   help, Moon_sim_image
  ;;   >  DOUBLE    = Array[800, 800, 160]

  ;;
  ;; QL file
  ;; 
  ofname_jpeg = outdir+'ql_'+file_basename(ofname_tiff,'.tiff') + '.jpg'
  tmp_image = out_hyper_image[*,*,40]
  byte_tmp_image = byte( tmp_image / (max(tmp_image) < 120) * 255d)
  write_jpeg,ofname_jpeg, byte_tmp_image

  return
end