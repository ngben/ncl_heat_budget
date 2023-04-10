load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; This script calculates the heat budget over the central equatorial
; IO. The formula for this is as in Cai et al. 2014 Extreme IODs

begin

  pi       = 4.0*atan(1.0)
  radians  = pi/180.0
  radius   = 6370949.0
  rr       = radius*radians

; array of models
  model = new(32,string)
  model(0) = "ACCESS1-0"
  model(1) = "ACCESS1-3"
  model(2) = "bcc-csm1-1"
  model(3) = "bcc-csm1-1-m"
  model(4) = "CanESM2"
  model(5) = "CCSM4"
  model(6) = "CESM1-BGC"
  model(7) = "CESM1-CAM5"
  model(8) = "CMCC-CESM"
  model(9) = "CMCC-CM"
  model(10) = "CMCC-CMS"
  model(11) = "CNRM-CM5"
  model(12) = "CSIRO-Mk3-6-0"
  model(13) = "FGOALS-g2"
  model(14) = "GFDL-CM3"
  model(15) = "GFDL-ESM2G"
  model(16) = "GFDL-ESM2M"
  model(17) = "GISS-E2-H"
  model(18) = "GISS-E2-R"
  model(19) = "HadGEM2-CC"
  model(20) = "HadGEM2-ES"
  model(21) = "IPSL-CM5A-LR"
  model(22) = "IPSL-CM5A-MR"
  model(23) = "IPSL-CM5B-LR"
  model(24) = "MIROC5"
  model(25) = "MIROC-ESM"
  model(26) = "MIROC-ESM-CHEM"
  model(27) = "MPI-ESM-LR"
  model(28) = "MPI-ESM-MR"
  model(29) = "MRI-CGCM3"
  model(30) = "NorESM1-ME"
  model(31) = "NorESM1-M"
  nmod = dimsizes(model)-1

;  period = "hist"
;  season = "SON"
;  restrict = "yes"
;  ySt = 191101
;  yEd = 200512
;  detrend = "yes"
  print("Enter period")
  period = systemfunc("read ANS; echo $ANS")
  print("Restrict period?")
  restrict = systemfunc("read ANS; echo $ANS")
  if (restrict.eq."yes") then
    print("Enter start year")
    ySt = stringtoint(systemfunc("read ANS; echo $ANS")+"01")
    print("Enter end year")
    yEd = stringtoint(systemfunc("read ANS; echo $ANS")+"12")
  end if
  print("Detrend the data? By default, data is NOT detrended")
  detrend = systemfunc("read ANS; echo $ANS")
if detrend.eq."no" then
  detrend_name = ""
else if detrend.eq."yes" then
  detrend_name = "_detrended"
end if
end if

do m = 0, nmod
  ful = systemfunc("ls /data/ng04l/cmip5_360x181/uo_"+model(m)+"_"+period+"*.nc")
  fvl = systemfunc("ls /data/ng04l/cmip5_360x181/vo_"+model(m)+"_"+period+"*.nc")
  fwl = systemfunc("ls /data/ng04l/cmip5_360x181/wo_"+model(m)+"_"+period+"*.nc")
  ftl = systemfunc("ls /data/ng04l/cmip5_360x181/thetao_"+model(m)+"_"+period+"*.nc")
  fswul = systemfunc(" ls /data/ng04l/cmip5_360x181/rsus_"+model(m)+"_"+period+"*.nc")
  fswdl = systemfunc(" ls /data/ng04l/cmip5_360x181/rsds_"+model(m)+"_"+period+"*.nc")
  flwul = systemfunc(" ls /data/ng04l/cmip5_360x181/rlus_"+model(m)+"_"+period+"*.nc")
  flwdl = systemfunc(" ls /data/ng04l/cmip5_360x181/rlds_"+model(m)+"_"+period+"*.nc")
  flhfl = systemfunc(" ls /data/ng04l/cmip5_360x181/hfls_"+model(m)+"_"+period+"*.nc")
  fshfl = systemfunc(" ls /data/ng04l/cmip5_360x181/hfss_"+model(m)+"_"+period+"*.nc")

  if ismissing(ful).or.ismissing(fvl).or.ismissing(fwl).or.ismissing(ftl)\
    .or.ismissing(fswul).or.ismissing(fswdl).or.ismissing(flwul).or.ismissing(flwdl)\
    .or.ismissing(flhfl).or.ismissing(fshfl) then
    continue
  end if
  fu = addfile(ful(0), "r")
  fv = addfile(fvl(0), "r")
  fw = addfile(fwl(0), "r")
  ft = addfile(ftl(0), "r")

if (restrict.eq."yes") then
  time_ = fu->time
  YYYY  = cd_calendar(time_,-1)
  delete(time_)
  tSt   = ind(YYYY.eq.ySt)
  tEd   = ind(YYYY.eq.yEd)
  uo = fu->uo(tSt:tEd,{0:100},:,:)
else if (restrict.eq."no") then
  uo = fu->uo(:,{0:100},:,:)
end if
end if
if detrend.eq."yes" then
  uod = dtrend_msg_n(uo&time, uo, False, False, 0) 
  copy_VarMeta(uo,uod)
  delete(uo)
  uo = uod
  delete(uod)
end if
  uob = clmMonTLLL(uo) ; b is for bar (i.e. the mean)
  uop = calcMonAnomTLLL(uo, uob) ; p is for prime (i.e. the anomaly)
  delete(uo)

; calculate mean and anomalous meridional current
if (restrict.eq."yes") then
  vo = fv->vo(tSt:tEd,{0:100},:,:)
else if (restrict.eq."no") then
  vo = fv->vo(:,{0:100},:,:)
end if
end if
if detrend.eq."yes" then
  vod = dtrend_msg_n(vo&time, vo, False, False, 0) 
  copy_VarMeta(vo,vod)
  delete(vo)
  vo = vod
  delete(vod)
end if
  vob = clmMonTLLL(vo)
  vop = calcMonAnomTLLL(vo, vob)
  delete(vo)

; calculate mean and anomalous vertical current, make sure sign is correct!
; wo data in /data/ng04l/cmip5_360x181 should be the correct sign
; (shoaling/upward current is positive!)
if (restrict.eq."yes") then
  wo = fw->wo(tSt:tEd,{0:100},:,:)
else if (restrict.eq."no") then
  wo = fw->wo(:,{0:100},:,:)
end if
end if
if detrend.eq."yes" then
  wod = dtrend_msg_n(wo&time, wo, False, False, 0) 
  copy_VarMeta(wo,wod)
  delete(wo)
  wo = wod
  delete(wod)
end if
  wob = clmMonTLLL(wo)
  wop = calcMonAnomTLLL(wo, wob)
  delete(wo)

; calculate mean temperature gradients (e.g. dTbdx)
if (restrict.eq."yes") then
  tmp = ft->thetao(tSt:tEd,{0:100},:,:)
else if (restrict.eq."no") then
  tmp = ft->thetao(:,{0:100},:,:)
end if
end if
if detrend.eq."yes" then
  tmpd = dtrend_msg_n(tmp&time, tmp, False, False, 0) 
  copy_VarMeta(tmp,tmpd)
  delete(tmp)
  tmp = tmpd
  delete(tmpd)
end if
  tmp_clm = clmMonTLLL(tmp)

  month  = tmp_clm&month
  lat    = tmp_clm&lat
  lon    = tmp_clm&lon
  lev    = tmp_clm&lev
  nlev   = dimsizes(lev)
  nlat   = dimsizes(lat)
  nlon   = dimsizes(lon)
  nmonth = dimsizes(month)
  dlat   = (max(lat)-min(lat))/(nlat-1) ; this is converted to metres later on
  dlon   = ((max(lon)-min(lon))/(nlon-1))*radians ; this is converted to metres later on

; calculate zonal mean temperature gradient
  dTbdx = new(dimsizes(tmp_clm),typeof(tmp_clm))
  copy_VarMeta(tmp_clm,dTbdx)
do i = 0, 11
  tmp_c = tmp_clm(i,:,:,:)
  dTbdx_= new(dimsizes(tmp_c),typeof(tmp_c),tmp_c@_FillValue)
  do nl = 0, nlat-1
    dx = doubletofloat(radius*cos(radians*lat(nl))*dlon) ;no need radians here as dlon is multiplied (by radians) already (see above), the other radians is for lat
    dTbdx_(:,nl:nl,:) = center_finite_diff_n(tmp_c(:,nl:nl,:), dx, False, 0, 2)
  end do
  copy_VarMeta(tmp_c, dTbdx_)
  dTbdx(i,:,:,:) = dTbdx_
  delete(dTbdx_)
  delete(tmp_c)
end do
  delete(i)
  dTbdx@units = "K/m"
  dTbdx@name  = "zonal mean temperature gradient"
  dTbdx@info  = "As the anomaly is calculated from the historical mean, the mean temperature used here is for the historical period"
  dTbdx@long_name = "zonal_mean_temperature_gradient"  
  print(model(m))

; calculate meridional mean temperature gradient
  dTbdy = new(dimsizes(tmp_clm),typeof(tmp_clm))
  copy_VarMeta(tmp_clm,dTbdy)
do i = 0, 11
  tmp_c = tmp_clm(i,:,:,:)
  dTbdy_= new(dimsizes(tmp_c), typeof(tmp_c), tmp_c@_FillValue)
  do ml=0,nlon-1
    dy = doubletofloat(dlat*radius*radians)
    dTbdy_(:,:,ml:ml) = center_finite_diff_n(tmp_c(:,:,ml:ml),dy,False,0,1)
  end do
  copy_VarMeta(tmp_c, dTbdy_)
  dTbdy(i,:,:,:) = dTbdy_
  delete(dTbdy_)
  delete(tmp_c)
end do
  delete(i)
  dTbdy@units = "K/m"
  dTbdy@name  = "meridional mean temperature gradient"
  dTbdy@info  = "As the anomaly is calculated from the historical mean, the mean temperature used here is for the historical period"
  dTbdy@long_name = "meridional_mean_temperature_gradient"  

; calculate vertical mean temperature gradient
;  dTbdz = new(dimsizes(tmp_clm),typeof(tmp_clm))
;  copy_VarMeta(tmp_clm,dTbdz)
;do i = 0, 11
;  tmp_c = tmp_clm(i,:,:,:)
;  dTbdz_= new(dimsizes(tmp_c), typeof(tmp_c), tmp_c@_FillValue)
;  z     = doubletofloat(lev)
;  dTbdz_= center_finite_diff_n(tmp_c, z, False, 0, 0)*(-1)
;  copy_VarMeta(tmp_c, dTbdz_)
;  dTbdz(i,:,:,:) = dTbdz_
;  delete(dTbdz_)
 ; delete(tmp_c)
;end do
;  delete(i)
;  dTbdz@units = "K/m"
;  dTbdz@name  = "vertical mean temperature gradient"
;  dTbdz@info  = "As the anomaly is calculated from the historical mean, the mean temperature used here is for the historical period"
;  dTbdz@long_name = "vertical_mean_temperature_gradient"  

; calculate vertical mean temperature gradient using Agus' method
  dTbdz = dim_avg_n_Wrap(tmp_clm(:,{0:50},:,:),1) ; avg. over top 50m
  T60_clm_av = tmp_clm(:,{60},:,:)
  dTbdz = (dTbdz - T60_clm_av)/(60-(50/2))
  delete(T60_clm_av)
  dTbdz@units = "K/m"
  dTbdz@name  = "vertical mean temperature gradient"
  dTbdz@info  = "dTdz is the average temperature in the upper 50m minus the temp at 60m depth. As the anomaly is calculated from the historical mean, the mean temperature used here is for the historical period"
  dTbdz@long_name = "vertical_mean_temperature_gradient"  

  delete(lat)
  delete(lon)
  delete(lev)
  delete(nlev)
  delete(nlat)
  delete(nlon)
  delete(dlat)
  delete(dlon)

; calculate anomalous temperature gradients (e.g. dTpdx)
  tmp_anom = calcMonAnomTLLL(tmp,tmp_clm)
  delete(tmp)
  delete(tmp_clm)
  time= tmp_anom&time
  lat = tmp_anom&lat
  lon = tmp_anom&lon
  lev = tmp_anom&lev
  nlev= dimsizes(lev)
  nlat= dimsizes(lat)
  nlon= dimsizes(lon)
  ntime= dimsizes(time)
  dlat = (max(lat)-min(lat))/(nlat-1)
  dlon = ((max(lon)-min(lon))/(nlon-1))*radians
  dtime= (max(time)-min(time))/(ntime-1)*86400 ; 86400 = number of seconds per day
; max(time) minus min(time) divided by ntime is the number of days 
; between each point (assuming it is evenly spaced).

; calcuate zonal anomalous temperature gradient
  print("Calculating zonal anomalous temperature gradient "+systemfunc("date"))
  dTpdx = new(dimsizes(tmp_anom),typeof(tmp_anom))
  copy_VarMeta(tmp_anom,dTpdx)
do i = 0, ntime-1
  tmp_a = tmp_anom(i,:,:,:)
  dTpdx_= new(dimsizes(tmp_a),typeof(tmp_a),tmp_a@_FillValue)
  do nl = 0, nlat-1
    dx = doubletofloat(radius*cos(radians*lat(nl))*dlon)
    dTpdx_(:,nl:nl,:) = center_finite_diff_n(tmp_a(:,nl:nl,:), dx, False, 0, 2)
  end do
  copy_VarMeta(tmp_a, dTpdx_)
  dTpdx(i,:,:,:) = dTpdx_
  delete(dTpdx_)
  delete(tmp_a)
end do
  delete(i)
  dTpdx@units = "K/m"
  dTpdx@name  = "zonal anomalous temperature gradient"
  dTpdx@info  = "The monthly anomaly for the RCP8.5 period is calculated from the historical mean"
  dTpdx@long_name = "zonal_anomalous_temperature_gradient"  

; calculate meridional anomalous temperature gradient
  print("Calculating meridional anomalous temperature gradient "+systemfunc("date"))
  dTpdy = new(dimsizes(tmp_anom),typeof(tmp_anom))
  copy_VarMeta(tmp_anom,dTpdy)
do i = 0, ntime-1
  tmp_a = tmp_anom(i,:,:,:)
  dTpdy_= new(dimsizes(tmp_a), typeof(tmp_a), tmp_a@_FillValue)
  do ml=0,nlon-1
    dy = doubletofloat(dlat*radius*radians)
    dTpdy_(:,:,ml:ml) = center_finite_diff_n(tmp_a(:,:,ml:ml),dy,False,0,1)
  end do
  copy_VarMeta(tmp_a, dTpdy_)
  dTpdy(i,:,:,:) = dTpdy_
  delete(dTpdy_)
  delete(tmp_a)
end do
  delete(i)
  dTpdy@units = "K/m"
  dTpdy@name  = "meridional anomalous temperature gradient"
  dTpdy@info  = "The monthly anomaly for the RCP8.5 period is calculated from the historical mean"
  dTpdy@long_name = "meridional_anomalous_temperature_gradient"

; calculate vertical anomalous temperature gradient, CHECK SIGN
  print("Calculating vertical anomalous temperature gradient "+systemfunc("date"))
;  dTpdz = new(dimsizes(tmp_anom),typeof(tmp_anom))  
;  copy_VarMeta(tmp_anom,dTpdz)
;do i = 0, ntime-1
;  tmp_a = tmp_anom(i,:,:,:)
;  dTpdz_= new(dimsizes(tmp_a), typeof(tmp_a), tmp_a@_FillValue)
;  z     = doubletofloat(lev)
;  dTpdz_= center_finite_diff_n(tmp_a, z, False, 0, 0)*(-1)
;  copy_VarMeta(tmp_a, dTpdz_)
;  dTpdz(i,:,:,:) = dTpdz_
;  delete(dTpdz_)
;  delete(tmp_a)
;end do
;  delete(i)
;  dTpdz@units = "K/m"
;  dTpdz@name  = "vertical anomalous temperature gradient"
;  dTpdz@info  = "The monthly anomaly for the RCP8.5 period is calculated from the historical mean"
;  dTpdz@long_name = "vertical_anomalous_temperature_gradient"

; calculate vertical anomalous temperature gradient using Agus' method
  dTpdz = dim_avg_n_Wrap(tmp_anom(:,{0:50},:,:),1) ; avg. over top 50m
  T60_anom_av= tmp_anom(:,{60},:,:)
  dTpdz = (dTpdz - T60_anom_av)/(60-(50/2))
  delete(T60_anom_av)
  dTpdz@units = "K/m"
  dTpdz@name  = "vertical anomalous temperature gradient"
  dTpdz@info  = "dTdz is the average temperature in the upper 50m minus the temp at 60m depth. As the anomaly is calculated from the historical mean, the mean temperature used here is for the historical period"
  dTpdz@long_name = "vertical_anomalous_temperature_gradient"  

; calculate anomalous temperature tendency (dT'/dt)
  print("Calculating anomalous temperature tendency "+systemfunc("date"))
  dTpdt = new(dimsizes(tmp_anom),typeof(tmp_anom))
  copy_VarMeta(tmp_anom,dTpdt)
  t     = doubletofloat(time)*86400 ; convert time from days since to seconds since.
; need a do loop? What do you loop over?
  dTpdt = center_finite_diff_n(tmp_anom, t, False, 0, 0) ; is this correct?
  delete(tmp_anom)
  dTpdt@units = "K/sec"
  dTpdt@name  = "anomalous temperature tendency"
  dTpdt@info  = "The monthly anomaly for the RCP8.5 period is calculated from the historical mean"
  dTpdt@long_name = "anomalous_temperature_tendency"

  print("Begin calculating heat fluxes "+systemfunc("date"))
; calculate net heat flux (Q) from air into sea i.e. downwards. Q = SW + LW + LH + SH
; swu = shortwave up, swd = shortwave down, lwu = longwave up, lwd = longwave down
; lhf = latent heat flux, shf = sensible heat flux. CHECK SIGN! Q should be positive down.
  fswu = addfile(fswul(0),"r")
  fswd = addfile(fswdl(0),"r")
  flwu = addfile(flwul(0),"r")
  flwd = addfile(flwdl(0),"r")
  flhf = addfile(flhfl(0),"r")
  fshf = addfile(fshfl(0),"r")
; calculate mean and anomalous downward shortwave radiation (sw)
if (restrict.eq."yes") then
  swu = fswu->rsus(tSt:tEd,:,:) ; rsus is upward positive
  swd = fswd->rsds(tSt:tEd,:,:) ; rsds is downward positive
else if (restrict.eq."no") then
  swu = fswu->rsus(:,:,:)
  swd = fswd->rsds(:,:,:)
end if
end if
  sw = swd
  sw = swd - swu ;this is rsds minus rsus
  delete(swd)
  delete(swu)
if detrend.eq."yes" then
  swd = dtrend_msg_n(sw&time, sw, False, False, 0) 
  copy_VarMeta(sw,swd)
  delete(sw)
  sw = swd
end if
  sw_clm      = clmMonTLL(sw) ; sw is downward positive
  sw_anom     = calcMonAnomTLL(sw, sw_clm)
  sw_anom@name = "net incoming shortwave radiation anomaly"

; calculate mean and anomalous downward longwave radiation (lw)
if (restrict.eq."yes") then
  lwu = flwu->rlus(tSt:tEd,:,:) ; rlus is upward positive
  lwd = flwd->rlds(tSt:tEd,:,:) ; rlds is downward positive
else if (restrict.eq."no") then
  lwu = flwu->rlus(:,:,:)
  lwd = flwd->rlds(:,:,:)
end if
end if
  lw = lwd
  lw = lwd - lwu ;this is rlds minus rlus
  delete(lwd)
  delete(lwu)
if detrend.eq."yes" then
  lwd = dtrend_msg_n(lw&time, lw, False, False, 0) 
  copy_VarMeta(lw,lwd)
  delete(lw)
  lw = lwd
  delete(lwd)
end if
  lw_clm      = clmMonTLL(lw) ; lw is downward positive
  lw_anom     = calcMonAnomTLL(lw, lw_clm)
  lw_anom@name = "net incoming longwave radiation anomaly"

; calculate mean and anomalous latent heat flux (lhf)
if (restrict.eq."yes") then
  lhf = flhf->hfls(tSt:tEd,:,:) ; hfls is upward positive
else if (restrict.eq."no") then
  lhf = flhf->hfls(:,:,:)
end if
end if
if detrend.eq."yes" then
  lhfd = dtrend_msg_n(lhf&time, lhf, False, False, 0) 
  copy_VarMeta(lhf,lhfd)
  delete(lhf)
  lhf = lhfd
  delete(lhfd)
end if
  lhf_clm      = clmMonTLL(lhf)
  lhf_anom     = calcMonAnomTLL(lhf, lhf_clm)

; calculate mean and anomalous sensible heat flux (shf)
if (restrict.eq."yes") then
  shf = fshf->hfss(tSt:tEd,:,:) ; hfss is upward positive
else if (restrict.eq."no") then
  shf = fshf->hfss(:,:,:)
end if
end if
if detrend.eq."yes" then
  shfd = dtrend_msg_n(shf&time, shf, False, False, 0) 
  copy_VarMeta(shf,shfd)
  delete(shf)
  shf = shfd
  delete(shfd)
end if
  shf_clm      = clmMonTLL(shf)
  shf_anom     = calcMonAnomTLL(shf, shf_clm)

; Qnet = SW - LW - SHF - LHF
; Qnet is positive downwards
; SW is positive downwards
; LW is positive UPWARDS
; SHF is positive UPWARDS
; LHF is positive UPWARDS
; OR Qnet = SW + LW - SHF - LHF if LW is positive downwards.

  Qnet_anom = sw_anom
  Qnet_anom@long_name = "Incoming Net Heat Flux anomaly"
  Qnet_anom@standard_name = "incoming_net_heat_flux_anomaly"
  Qnet_anom = sw_anom + lw_anom - lhf_anom - shf_anom

; divide Q by density of water, specific heat of water and the depth of the mixed layer (50m)
  Qnet_anom = Qnet_anom/(1026*3986*50) ; this should give units of degrees Kelvin/second
; 3996 instead of 3986 may be a more suitable specific heat capacity of water
  Qnet_anom@units = "K/m"
  Qnet_anom@info = "Net heat flux has been divided by the density of water, the heat "

; NOTE: cannot calculate residual without averaging over depth. As such,
; calculation of the residual is done when examining the heat budget over a
; specific area (e.g. the CEIO) and is not output in this script.

  delete(time)
  delete(ntime)
  delete(lat)
  delete(lon)
  delete(lev)
  time= uop&time
  lat = uop&lat
  lon = uop&lon
  lev = uop&lev
  ntime= dimsizes(time)
  nlev= dimsizes(lev)
  nlat= dimsizes(lat)
  nlon= dimsizes(lon)

; output data
  print("Outputting results "+systemfunc("date"))
  filo = model(m)+"_"+period+"_heat_budget"+detrend_name+".nc"
  lso = systemfunc("ls -d ../analysis_data/heat_budget"+detrend_name)
  if ismissing(lso) then
    mko = systemfunc("mkdir ../analysis_data/heat_budget"+detrend_name)
  end if
  system("mv ../analysis_data/heat_budget"+detrend_name+"/"+filo+" ../analysis_data/heat_budget"+detrend_name+"/"+filo+"_tmp")
  fout  = addfile("../analysis_data/heat_budget"+detrend_name+"/"+filo,"c")
  
  setfileoption(fout,"DefineMode",True)
  fAtt       = True
  fAtt@title = model(m)+" global heat budget analysis"
  fAtt@source_files = ful+"\n"+fvl+"\n"+fwl+"\n"+ftl+"\n"+fswul+"\n"+fswdl+"\n"+flwul\
  +"\n"+flwdl+"\n"+flhfl+"\n"+fshfl
  fAtt@Conventions  = "None"
  fAtt@creation_date= systemfunc("date")
if detrend.eq."no" then
  fAtt@info = "Data has not been detrended"
else if detrend.eq."yes" then
  fAtt@info = "Data has been detrended"
end if
end if
  fileattdef(fout,fAtt)
  
  dimNames = (/"time","month","lev","lat","lon"/)
  dimSizes = (/-1,nmonth,nlev,nlat,nlon/)
  dimUnlim = (/True,False,False,False,False/)
  filedimdef(fout,dimNames,dimSizes,dimUnlim)
  
  filevardef(fout,"time",typeof(time),getvardims(time))
  filevardef(fout,"month",typeof(month),getvardims(month)) 
  filevardef(fout, "lev" ,typeof(lev),getvardims(lev))
  filevardef(fout, "lat" ,typeof(lat),getvardims(lat))    
  filevardef(fout, "lon" ,typeof(lon),getvardims(lon))
  filevardef(fout,"uop",typeof(uop),getvardims(uop))
  filevardef(fout,"uob",typeof(uob),getvardims(uob))
  filevardef(fout,"vop",typeof(vop),getvardims(vop))
  filevardef(fout,"vob",typeof(vob),getvardims(vob))
  filevardef(fout,"wop",typeof(wop),getvardims(wop))
  filevardef(fout,"wob",typeof(wob),getvardims(wob))
  filevardef(fout,"dTpdx",typeof(dTpdx),getvardims(dTpdx))
  filevardef(fout,"dTbdx",typeof(dTbdx),getvardims(dTbdx))
  filevardef(fout,"dTpdy",typeof(dTpdy),getvardims(dTpdy))
  filevardef(fout,"dTbdy",typeof(dTbdy),getvardims(dTbdy))
  filevardef(fout,"dTpdz",typeof(dTpdz),getvardims(dTpdz))
  filevardef(fout,"dTbdz",typeof(dTbdz),getvardims(dTbdz))
;  filevardef(fout,"dTpdz_agus",typeof(dTpdz_agus),getvardims(dTpdz_agus))
;  filevardef(fout,"dTbdz_agus",typeof(dTbdz_agus),getvardims(dTbdz_agus))
  filevardef(fout,"dTpdt",typeof(dTpdt),getvardims(dTpdt))
  filevardef(fout,"Qnet_anom",typeof(Qnet_anom),getvardims(Qnet_anom))
  filevardef(fout,"sw_anom",typeof(sw_anom),getvardims(sw_anom))
  filevardef(fout,"lw_anom",typeof(lw_anom),getvardims(lw_anom))
  filevardef(fout,"shf_anom",typeof(shf_anom),getvardims(shf_anom))
  filevardef(fout,"lhf_anom",typeof(lhf_anom),getvardims(lhf_anom))

  filevarattdef(fout,"time",time)
  filevarattdef(fout,"month",month)
  filevarattdef(fout,"lev"  ,lev)                     ; copy lev 
  filevarattdef(fout,"lat"  ,lat)                     ; copy lat attributes
  filevarattdef(fout,"lon"  ,lon)                     ; copy lon attributes
  filevarattdef(fout,"uop",uop)
  filevarattdef(fout,"uob",uob)
  filevarattdef(fout,"vop",vop)
  filevarattdef(fout,"vob",vob)
  filevarattdef(fout,"wop",wop)
  filevarattdef(fout,"wob",wob)
  filevarattdef(fout,"dTpdx",dTpdx)
  filevarattdef(fout,"dTbdx",dTbdx)
  filevarattdef(fout,"dTpdy",dTpdy)
  filevarattdef(fout,"dTbdy",dTbdy)
  filevarattdef(fout,"dTpdz",dTpdz)
  filevarattdef(fout,"dTbdz",dTbdz)
;  filevarattdef(fout,"dTpdz_agus",dTpdz_agus)
;  filevarattdef(fout,"dTbdz_agus",dTbdz_agus)
  filevarattdef(fout,"dTpdt",dTpdt)
  filevarattdef(fout,"Qnet_anom",Qnet_anom)
  filevarattdef(fout,"sw_anom",sw_anom)
  filevarattdef(fout,"lw_anom",lw_anom)
  filevarattdef(fout,"shf_anom",shf_anom)
  filevarattdef(fout,"lhf_anom",lhf_anom)

  setfileoption(fout,"DefineMode",False)
  fout->uop   = (/uop/)
  fout->uob   = (/uob/)
  fout->vop   = (/vop/)
  fout->vob   = (/vob/)
  fout->wop   = (/wop/)
  fout->wob   = (/wob/)
  fout->dTpdx = (/dTpdx/)
  fout->dTpdy = (/dTpdy/)
  fout->dTpdz = (/dTpdz/)
  fout->dTbdx = (/dTbdx/)
  fout->dTbdy = (/dTbdy/)
  fout->dTbdz = (/dTbdz/)
  fout->dTpdt = (/dTpdt/)
;  fout->dTbdz_agus = (/dTbdz_agus/)
;  fout->dTpdz_agus = (/dTpdz_agus/)
  fout->Qnet_anom  = (/Qnet_anom/)
  fout->sw_anom    = (/sw_anom/)
  fout->lw_anom    = (/lw_anom/)
  fout->shf_anom   = (/shf_anom/)
  fout->lhf_anom   = (/lhf_anom/)
  fout->time       = (/time/)
  fout->month      = (/month/)
  fout->lev        = (/lev/)    
  fout->lat        = (/lat/)
  fout->lon        = (/lon/)
  curdate         = systemfunc("date +%Y.%m.%d")

  delete(time)
  delete(uop)
  delete(uob)
  delete(vop)
  delete(vob)
  delete(wop)
  delete(wob)
  delete(dTpdx)
  delete(dTbdx)
  delete(dTpdy)
  delete(dTbdy)
  delete(dTpdz)
  delete(dTbdz)
;  delete(dTpdz_agus)
;  delete(dTbdz_agus)
  delete(dTpdt)
  delete(Qnet_anom)
  delete(sw_anom)
  delete(lw_anom)
  delete(shf_anom)
  delete(lhf_anom)
  delete(lon)
  delete(lat)
  delete(lev)
if restrict.eq."yes" then
  delete(YYYY)
  delete(tSt)
  delete(tEd)
end if
  delete(t) ;used to calculate dTpdt
  delete(month)

end do
end
