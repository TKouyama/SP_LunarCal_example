;; SP model data

function channel_wavelength,arr_size=arr_size

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
    1620.099976, 1628.099976, 1636.099976, 1644.199951, $
    ;; N2 ;;
    1709.8, 1717.6, 1725.6, 1733.7, 1742.0, 1749.7, 1757.7, 1766.3, 1773.6, 1782.2, $
    1789.8, 1797.6, 1805.8, 1813.7, 1822.0, 1830.0, 1837.6, 1845.6, 1853.7, 1861.8, $
    1870.1, 1877.3, 1885.7, 1893.7, 1901.5, 1910.0, 1918.0, 1925.3, 1934.3, 1942.0, $
    1948.8, 1957.6, 1965.9, 1973.3, 1981.3, 1989.4, 1997.7, 2005.8, 2013.0, 2021.5, $
    2029.3, 2037.4, 2045.8, 2053.3 $
  ]
  
  c_size=size(wavelength)
  arr_size=c_size[1]
  return, wavelength
end
