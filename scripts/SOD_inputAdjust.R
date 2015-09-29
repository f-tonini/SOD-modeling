####### RASTER INPUT ADJUSTMENT ##############
#adjust for data errors in case it's over live trees!
host_sum_rast <- umca_rast + sese_rast + acma_rast + 
  aeca_rast + arme_rast + psme_rast + 
  lide_rast + oaks_rast

#if host_sum > live trees layer, THEN ADJUST SO THAT: host_sum = live trees
#the adjustment is done by redistributing densities based on the relative abundance
if(any(host_sum_rast[] > lvtree_rast[])){
  
  cond_true <- host_sum_rast[] > lvtree_rast[]
  umca_rast[cond_true] <- floor(lvtree_rast[cond_true] * (umca_rast[cond_true] / host_sum_rast[cond_true]))
  sese_rast[cond_true] <- floor(lvtree_rast[cond_true] * (sese_rast[cond_true] / host_sum_rast[cond_true]))
  acma_rast[cond_true] <- floor(lvtree_rast[cond_true] * (acma_rast[cond_true] / host_sum_rast[cond_true]))
  aeca_rast[cond_true] <- floor(lvtree_rast[cond_true] * (aeca_rast[cond_true] / host_sum_rast[cond_true]))
  arme_rast[cond_true] <- floor(lvtree_rast[cond_true] * (arme_rast[cond_true] / host_sum_rast[cond_true]))
  psme_rast[cond_true] <- floor(lvtree_rast[cond_true] * (psme_rast[cond_true] / host_sum_rast[cond_true]))
  lide_rast[cond_true] <- floor(lvtree_rast[cond_true] * (lide_rast[cond_true] / host_sum_rast[cond_true]))
  oaks_rast[cond_true] <- floor(lvtree_rast[cond_true] * (oaks_rast[cond_true] / host_sum_rast[cond_true]))
  writeRaster(umca_rast, "./layers/adjusted/UMCA_den_100m.img", format='HFA', datatype='INT1U', overwrite=TRUE)
  writeRaster(sese_rast, "./layers/adjusted/SESE3_den_100m.img", format='HFA', datatype='INT1U', overwrite=TRUE) 
  writeRaster(acma_rast, "./layers/adjusted/ACMA3_den_100m.img", format='HFA', datatype='INT1U', overwrite=TRUE) 
  writeRaster(aeca_rast, "./layers/adjusted/AECA_den_100m.img", format='HFA', datatype='INT1U', overwrite=TRUE) 
  writeRaster(arme_rast, "./layers/adjusted/ARME_den_100m.img", format='HFA', datatype='INT1U', overwrite=TRUE) 
  writeRaster(psme_rast, "./layers/adjusted/PSME_den_100m.img", format='HFA', datatype='INT1U', overwrite=TRUE) 
  writeRaster(lide_rast, "./layers/adjusted/LIDE3_den_100m.img", format='HFA', datatype='INT1U', overwrite=TRUE) 
  writeRaster(oaks_rast, "./layers/adjusted/OAKS_den_100m.img", format='HFA', datatype='INT1U', overwrite=TRUE) 
}
