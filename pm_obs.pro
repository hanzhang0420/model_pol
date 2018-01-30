pro pm_obs,incl,output,R_val,file_b_field

  ;0.01-1.0 dis
  kappa_p=double(4080.5)   ;at wavelength 10 mum   (v55(j)/m)*!pi*(d)^2
  kappa_l=double(2833.75)
  print,'incl of the disk',incl
  print,'R alignment eff',R_val
  wave=10.3*(10.0^4.0)
  file='/Users/hanzhang/polmodel/from_radmc/dt_cart_ab_aur.xdr'
  restore,filename=file_b_field;'/scratch/hanzhang/plot_model/b_field/b_com.xdr'

  nout=51
  pm_integrate,kappa_p,kappa_l,file,ms,incl,I_obs,U_obs,Q_obs,wave,nout,R_val
  print,'size of final image',nout
  print,'finish integration, saving file...'
  save,I_obs,U_obs,Q_obs,filename=output
end