pro angles_to_xyz,r,phi,theta,x,y,z
  DRADEG = 180.d0/!dpi
  x = r * cos(phi / DRADEG) * sin(theta / DRADEG)
  y = r * sin(phi / DRADEG) * sin(theta / DRADEG)
  z = r * cos(theta / DRADEG)
  return
end

pro pm_coord_regrid,filein_grid,filein_dt,fileout,ncube,bl,bu,LINEAR=linear

  ; actually rc tc pc in the problem_setup_yso.pro
  ; regrid from spherical coordinate to cartesian coordinate
  ; input:
  ;      filein_grid: grid output file  from radmc
  ;      filein_dt: density and tempertaure output file from radmc
  ;      fileout: file to be saved (dens, temp)
  ;      ngrid: number of grids in the cube
  ;      boundary: bl(lower) bu(upper) 
  
  restore,filename=filein_grid
  
  nr=n_elements(grid_r_s)
  ntheta=n_elements(grid_theta_s)
  nphi=n_elements(grid_phi_s)

  print,'size of coordinate',nr,ntheta,nphi
  au=1.49597871*10.0^13.0
  grid_r=grid_r_s/au ;change the unit to au

  x_t=fltarr(nr,ntheta,nphi)
  y_t=fltarr(nr,ntheta,nphi)
  z_t=fltarr(nr,ntheta,nphi)

  ;from spherical to cartesian
  for i=0,nr-1 do begin
    for j=0,ntheta-1 do begin
      for k=0,nphi-1 do begin
        angles_to_xyz,grid_r[i],grid_phi_s[k]*180/!pi,grid_theta_s[j]*180/!pi,x,y,z
        x_t[i,j,k]=x
        y_t[i,j,k]=y
        z_t[i,j,k]=z

      endfor
    endfor
  endfor  
   
  restore,filename=filein_dt
  ;format of dens_temp_s_1: dens, temp of the three dimensional structure
  ;!important reform to one-dimensional
  x_t_r=reform(x_t,nr*ntheta*nphi)
  y_t_r=reform(y_t,nr*ntheta*nphi)
  z_t_r=reform(z_t,nr*ntheta*nphi)
  dens_r=reform(dens_s,nr*ntheta*nphi)
  temp_r=reform(temp_s,nr*ntheta*nphi)

  print,minmax(dens_r),minmax(temp_r)

  n_grid=n_elements(x_t_r)
  print,'number of x',n_grid
  ;print,minmax(dens_r),minmax(temp_r)
 
;  z_c=[-rotate(arrgen(bl,bu,/log,nstep=ncube/2.0),2),arrgen(bl,bu,/log,nstep=ncube/2.0)]
;  x_c=[-rotate(arrgen(bl,bu,/log,nstep=ncube/2.0),2),arrgen(bl,bu,/log,nstep=ncube/2.0)]
;  y_c=[-rotate(arrgen(bl,bu,/log,nstep=ncube/2.0),2),arrgen(bl,bu,/log,nstep=ncube/2.0)]
  
  if keyword_set(LINEAR) then begin
    print,'linear spaced array'
    z_c=findgen(ncube+1)*(2.0*bu)/(ncube+1.0-1.0)-bu
    x_c=findgen(ncube+1)*(2.0*bu)/(ncube+1.0-1.0)-bu
    y_c=findgen(ncube+1)*(2.0*bu)/(ncube+1.0-1.0)-bu
  endif

  x_new=(x_c[0:ncube-1]+x_c[1:ncube])/2.0
  y_new=(y_c[0:ncube-1]+y_c[1:ncube])/2.0
  z_new=(z_c[0:ncube-1]+z_c[1:ncube])/2.0

  n_new_grid=n_elements(z_new)
  
  print,'new_grid',n_new_grid,n_new_grid,n_new_grid
  
  regrid_dens=fltarr(n_new_grid,n_new_grid,n_new_grid)
  regrid_temp=fltarr(n_new_grid,n_new_grid,n_new_grid)


  for j=0,n_new_grid-1 do begin
    
    count=0l
    x_s=fltarr(1000000) & y_s=fltarr(1000000)
    dens_s=fltarr(1000000) & temp_s=fltarr(1000000)

    diff=z_t_r-z_new[j]
    

for i=0,n_grid-1 do begin
  if abs(z_t_r[i]-z_new[j]) le 0.2  then begin

    x_s[count]=x_t_r[i]
    y_s[count]=y_t_r[i]
    dens_s[count]=dens_r[i]
    temp_s[count]=temp_r[i]
    count=count+1
  endif

endfor
print,count,z_new[j]

x=x_s[0:count-1]
y=y_s[0:count-1]
print,minmax(abs(x)),minmax(abs(y))
dens=dens_s[0:count-1]
temp=temp_s[0:count-1]
TRIANGULATE, x, y, tr, b
res_dens=TRIGRID(x, y, dens, tr,xout=x_new,yout=y_new)
res_temp=TRIGRID(x, y, temp, tr,xout=x_new,yout=y_new)
print,minmax(res_dens),minmax(dens),minmax(res_temp),minmax(temp)
;  stop
regrid_dens[*,*,j]=res_dens
regrid_temp[*,*,j]=res_temp
   
    endfor
  
  save,x_c,x_new,regrid_dens,regrid_temp,filename=fileout
 
end

