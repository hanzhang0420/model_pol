pro pm_ms,para1,para2,bz_c,b_x,b_y,b_z,ncube

  z_c=findgen(ncube+1)*600.0/(ncube-1)-300.0
  x_c=findgen(ncube+1)*600.0/(ncube-1)-300.0
  y_c=findgen(ncube+1)*600.0/ncube-1-300.0

  x_new=(x_c[0:ncube-1]+x_c[1:ncube])/2.0
  y_new=(y_c[0:ncube-1]+y_c[1:ncube])/2.0
  z_new=(z_c[0:ncube-1]+z_c[1:ncube])/2.0
  Undefine,z_c,x_c,y_c
  n_grid=n_elements(z_new)
  print,n_grid
  b_x=fltarr(n_grid,n_grid,n_grid)
  b_y=fltarr(n_grid,n_grid,n_grid)
  b_z=fltarr(n_grid,n_grid,n_grid)


  for nx=0,n_grid-1 do begin
    for ny=0,n_grid-1 do begin
      for nz=0,n_grid-1 do begin

        r=sqrt(x_new[nx]^2+y_new[ny]^2)
        theta=atan(y_new[ny],x_new[nx])

        if r EQ 0 then begin
          B_r=0 & B_phi=0
        endif else begin
          ;  r=sqrt(x_new[nx]^2+y_new[ny]^2)
          ;   theta=atan(y_new[ny],x_new[nx])
          B_z[nx,ny,nz]=bz_c
          B_r=para1*(z_new[nz]/r)*B_z[nx,ny,nz]
          B_phi=para2*(z_new[nz]/r)*B_z[nx,ny,nz]

        endelse
        b_y[nx,ny,nz]=B_r*cos(theta)-B_phi*sin(theta)
        b_x[nx,ny,nz]=B_r*sin(theta)+B_phi*cos(theta)
      endfor
    endfor
  endfor

 ms=CREATE_STRUCT('b_x',b_x,'b_y',b_y,'b_z',b_z)
 return
end

pro write_ms
;;modelI in aitken
pm_ms,0.0,0.0,10.0,b_x,b_y,b_z
save,ms,filename='/Users/hanzhang/polmodel/b_field/b_I.xdr'
;
;;modelII
;b_field_orien_2,10.0,0,b_x,b_y,b_z
;save,b_x,b_y,b_z,filename='/scratch/hanzhang/plot_model/b_field/b_II.xdr'

end