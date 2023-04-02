# SP_LunarCal_example
Example code to read SP map model.

- example_of_read_sp_model.pro

Evaluating Moon surface radiance at a specified observation geometry based on map version of SP model.
The code will output radiance map (W/m2/sr/um) with tiff format as an example, and quick look of maps of incident, emission, and phase angles and some other information.

--- Coded by T. Kouyama (t.kouyama@aist.go.jp).

Please replace below path in the code with your directory name where SP files are contained.
L36. ifldname ='./'

All required functions are contained in the "pro" file.
(But if you find missed function(s), please inform me)

If you use this code for your research purpose,
please cite below papers in your documents and/or papers:

- Yokota et al., 2011: Lunar photometric properties at wavelengths 0.5-1.6μm acquired by SELENE Spectral Profiler and their dependency on local albedo and latitudinal zones, Icarus, 215, 639-660

- Ogohara et al., 2012: Automated cloud tracking system for the Akatsuki Venus Climate Orbiter data, Icarus, 217, 661-668

- Kouyama et al., 2016: Development of an application scheme for the SELENE/SP lunar reflectance model for radiometric calibration of hyperspectral and multispectral sensors, Planet. Space Sci., 124, 76-83

---
Structure of parameter files:

- 'avg_cube_1000s-7000s_selected_ip110225.img' :: Binary data, Double-precision float
SP model cube data (lon, lat, wavelength)
0.5x0.5 degree grid interval in longitude and latitude direction, and there are 160 channels in wavelength.
i.e. (lon, lat, wavelength) = (720, 360, 160)

  ;; Note: Definition of SP map coordinate  (Lon, Lat)
  ;; (0.25, 89.75), (0.75, 89.75), ... (359.75,89.75) ;; North pole
  ;; (0.25, 89.25), (0.75, 89.25), ... (359.75,89.25)
  ;; ...
  ;; (0.25, 0.25), (0.75, 0.25), ... (359.75,0.25)
  ;; (0.25, -0.25), (0.75, -0.25), ... (359.75,-0.25)
  ;; ...
  ;; (0.25, -89.25), (0.75, -89.25), ... (359.75,-89.25)
  ;; (0.25, -89.75), (0.75, -89.75), ... (359.75,-89.75) ;; South pole

- 'albedo_group_05x05.dat' :: Binary data, 2byte Integer
 (lon, lat) = (720, 360)
Albedo group map (this file is made by Kouyama based on Yokota's definition in Yokota et al 2011)
Note; Latitude direction may be opposit to the SP model definition, sorry..


- 'High_albedo_sel.txt'
- 'Mid_albedo_sel.txt'
- 'Low_albedo_sel.txt'
ASCII files. Parameters of i, e, alpha dependence for each albedo group.


- 'Gueymard.txt' :: ASCII file
Solar irradiance data from Gueymard model which is used in generating SP model


---

Toru Kouyama, 2019.04.29, modified 2020.10.12, 11.16
contact: t.kouyama@aist.go.jp
