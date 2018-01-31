pro res_yso

  file_b_field='/Users/hanzhang/polmodel/b_field/b_III_yso.xdr'
  output='/Users/hanzhang/polmodel/res_yso/b_III'
  incl=[0,30,60,80]
  for i=0,3 do begin
    pm_yso,incl[i],output+'_'+strtrim(incl[i],1)+'.xdr',0.25,file_b_field
  endfor

  print,'done yso'
  
end

pro res_yso_python


  filename='/Users/hanzhang/polmodel/res_yso/b_III'
  output='/Users/hanzhang/polmodel/res_yso/res_b_III'
    incl=[0,30,60,80]
  for i=0,3 do begin
    idl_stokes_incl,filename+'_'+strtrim(incl[i],1)+'.xdr',output+'_'+strtrim(incl[i],1)+'.xdr',51
  endfor
  
end