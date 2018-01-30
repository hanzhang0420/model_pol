function rot_incl_image,array,incl
  ;Angle of rotation in degrees clockwise.
  n_grid=size(array)
  ;rotate around x axis
  image_rot=fltarr(n_grid[1],n_grid[2],n_grid[3])
  for i=0,n_grid[1]-1 do begin
    origi=reform(array[i,*,*])
    image_rot[i,*,*]=rot(origi,incl,/INTERP)
  endfor
  return,image_rot
end

pro ms_b,b_x,b_y,b_z,incl,xi,gamma
  b_y_new=b_y*cos(incl*!pi/180)-b_z*sin(incl*!pi/180)
  b_z_new=b_z*cos(incl*!pi/180)+b_y*sin(incl*!pi/180)
  undefine,b_y,b_z
  xi=atan(b_y_new,b_x)-!pi/2.0
  gamma=atan(b_z_new,sqrt(b_y_new^2+b_x^2))
  undefine,b_x,b_y_new,b_z_new 
end


pro pm_integrate,kp,kl,file,ms,incl,I_obs,U_obs,Q_obs,wave,nout,R_val
  
  ;AU  = 1.49598d13     ; Astronomical Unit       [cm]
   ;at wavelength 10 mum

  restore,filename=file
  ; rot around x axis
  dens=rot_incl_image(regrid_dens,incl)
  Undefine,regrid_dens
  temp=rot_incl_image(regrid_temp,incl)
  Undefine,regrid_temp ;free memory

  n_grid=size(dens)
  print,'size of grid',n_grid[1],n_grid[2],n_grid[3]  
 
  ncube=n_grid[1]
  
  ms_b,ms.b_x,ms.b_y,ms.b_z,incl,xi,gam
  xi=xi+!pi/2.0
  
  R=double(R_val)
  
  I_i=dblarr(ncube,ncube)
  U_i=dblarr(ncube,ncube)
  Q_i=dblarr(ncube,ncube)
  

 S = CALL_EXTERNAL('pol_model.so', 'cal_integrate_', kp,kl,double(wave),double(dens),double(temp),R,$
  I_i,Q_i,U_i,fix(ncube),double(x_c),double(gam),double(xi))

  I_obs=rebin(I_i,nout,nout)
  Q_obs=rebin(Q_i,nout,nout)
  U_obs=rebin(U_i,nout,nout)
 print,minmax(I_obs)

end
