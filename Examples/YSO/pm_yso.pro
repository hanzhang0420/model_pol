pro pm_integrate_yso,kp,kl,file,ms,incl,I_obs,U_obs,Q_obs,wave,nout,R_val,x_c

  ;AU  = 1.49598d13     ; Astronomical Unit       [cm]
  ;at wavelength 10 mum

  restore,filename=file
  ; rot around x axis
  dens=rot_incl_image(dens,incl)
  Undefine,regrid_dens
  temp=rot_incl_image(temp,incl)
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


pro pm_yso,incl,output,R_val,file_b_field

;incl=0 & R_val=0.05
file='/Users/hanzhang/polmodel/from_radmc/dens_temp_yso.xdr'
;file_b_field='/Users/hanzhang/polmodel/b_field/b_III_yso.xdr'
;restore,'/Users/hanzhang/polmodel/from_radmc/grid_c.xdr'
;restore,'/Users/hanzhang/polmodel/from_radmc/dens_temp_yso.xdr'


;kappa_p=double(4080.5)   ;at wavelength 10 mum   (v55(j)/m)*!pi*(d)^2
;kappa_l=double(2833.75)

;at wavelength 10 mum 0.01-0.25
kappa_p=double(4083.2)   
kappa_l=double(2843.6)

print,'incl of the disk',incl
print,'R alignment eff',R_val
wave=10.3*(10.0^4.0)

restore,filename=file_b_field;'/scratch/hanzhang/plot_model/b_field/b_com.xdr'

ncube=153
x_c=findgen(ncube+1)*20000.0/(ncube)-10000.0

nout=51
pm_integrate_yso,kappa_p,kappa_l,file,ms,incl,I_obs,U_obs,Q_obs,wave,nout,R_val,x_c
print,'size of final image',nout
print,'finish integration, saving file...'
save,I_obs,U_obs,Q_obs,filename=output


end

