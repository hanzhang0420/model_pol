pro idl_stokes_incl,filename,output,ni
restore,filename=filename

;where does 3.9*10^9 come from???
I = rebin(I_obs*3.6146e11,ni,ni)/!pi;filter_image( rebin(I_obs,51,51), FWHM=[3.75], /ALL );
Q = rebin(Q_obs*3.6146e11,ni,ni)/!pi; filter_image( rebin(Q_obs,51,51), FWHM=[3.75], /ALL );
U = rebin(U_obs*3.6146e11,ni,ni)/!pi

;data_em=readfits('/scratch/hanzhang/ppdisk_ab_aur/image_ab_aur_1.2_em.fits',h)
;I=data_em
;Q=Q_obs*data_em/(I_obs+10.0^(-30))
;U=U_obs*data_em/(I_obs+10.0^(-30))

;data_em=fltarr(102,102,3)
;data_em[*,*,0]=I
;data_em[*,*,1]=Q
;data_em[*,*,2]=U

;I=rebin(I,51,51)
;Q=rebin(Q,51,51)
;U=rebin(U,51,51)

;I=filter_image( rebin(I,51,51), FWHM=[3.571], /ALL );
;Q=filter_image(rebin(Q,51,51), FWHM=[3.571], /ALL );
;U=filter_image(rebin(U,51,51), FWHM=[3.571], /ALL );

;atv,alog10(I)
help,I
n=size(I)
du=fltarr(n(1),n(2))
dq=fltarr(n(1),n(2))

m=5

tmpq_bin=(P_BIN_SIMPLE(Q,m))
tmpu_bin=(P_BIN_SIMPLE(U,m))
tmpu_i=(P_BIN_SIMPLE(I,m))

help,tmpu_i
;res=p_percent_theta(tmpu_bin/tmpu_i,tmpq_bin/tmpu_i,du,dq)
;p=percent_check(res.percent,1.0)
;theta=res.theta+90 ; by definition, the angle is east from north but for python to plot the image


ip=sqrt(Q^2+U^2)
theta=0.5*atan(tmpu_bin,tmpq_bin)+!pi/2.0
p=sqrt(tmpu_bin^2+tmpq_bin^2)/tmpu_i

theta_b=theta+!pi/2.0

save,p,theta,i,ip,filename=output


end






pro res_python


  filename='/Users/hanzhang/polmodel/res_image/b_I'
  output='/Users/hanzhang/polmodel/res_image/res_b_I'
  incl=[0,30,60,80]
  for i=0,3 do begin
    idl_stokes_incl,filename+'_'+strtrim(incl[i],1)+'.xdr',output+'_'+strtrim(incl[i],1)+'.xdr'
  endfor

  filename='/Users/hanzhang/polmodel/res_image/b_II'
  output='/Users/hanzhang/polmodel/res_image/res_b_II'
  incl=[0,30,60,80]
  for i=0,3 do begin
    idl_stokes_incl,filename+'_'+strtrim(incl[i],1)+'.xdr',output+'_'+strtrim(incl[i],1)+'.xdr'
  endfor
  
  filename='/Users/hanzhang/polmodel/res_image/b_III'
  output='/Users/hanzhang/polmodel/res_image/res_b_III'
  incl=[0,30,60,80]
  for i=0,3 do begin
    idl_stokes_incl,filename+'_'+strtrim(incl[i],1)+'.xdr',output+'_'+strtrim(incl[i],1)+'.xdr'
  endfor
;  filename='/scratch/hanzhang/plot_model/res_image/b_II'
;  output='/scratch/hanzhang/plot_model/res_image/res_b_II'
;  incl=[0,30,60,80]
;  for i=0,3 do begin
;    idl_stokes_incl,filename+'_'+strtrim(incl[i],1)+'.xdr',output+'_'+strtrim(incl[i],1)+'.xdr'
;  endfor
;
;
;  filename='/scratch/hanzhang/plot_model/res_image/b_III'
;  output='/scratch/hanzhang/plot_model/res_image/res_b_III'
;  incl=[0,30,60,80]
;  for i=0,3 do begin
;    idl_stokes_incl,filename+'_'+strtrim(incl[i],1)+'.xdr',output+'_'+strtrim(incl[i],1)+'.xdr'
;  endfor
;
;  filename='/scratch/hanzhang/plot_model/res_image/b_IV'
;  output='/scratch/hanzhang/plot_model/res_image/res_b_IV'
;  incl=[0,30,60,80]
;  for i=0,3 do begin
;    idl_stokes_incl,filename+'_'+strtrim(incl[i],1)+'.xdr',output+'_'+strtrim(incl[i],1)+'.xdr'
;  endfor
;
;
;  filename='/scratch/hanzhang/plot_model/res_image/b_V'
;  output='/scratch/hanzhang/plot_model/res_image/res_b_V'
;  incl=[0,30,60,80]
;  for i=0,3 do begin
;    idl_stokes_incl,filename+'_'+strtrim(incl[i],1)+'.xdr',output+'_'+strtrim(incl[i],1)+'.xdr'
;  endfor
;
;  filename='/scratch/hanzhang/plot_model/res_image/b_VI'
;  output='/scratch/hanzhang/plot_model/res_image/res_b_VI'
;  incl=[0,30,60,80]
;  for i=0,3 do begin
;    idl_stokes_incl,filename+'_'+strtrim(incl[i],1)+'.xdr',output+'_'+strtrim(incl[i],1)+'.xdr'
;  endfor
;
;  filename='/scratch/hanzhang/plot_model/res_image/b_VII'
;  output='/scratch/hanzhang/plot_model/res_image/res_b_VII'
;  incl=[0,30,60,80]
;  for i=0,3 do begin
;    idl_stokes_incl,filename+'_'+strtrim(incl[i],1)+'.xdr',output+'_'+strtrim(incl[i],1)+'.xdr'
;  endfor

end