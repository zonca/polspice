function read_cov, file, flatten=flatten

spawn,/sh,'wc -l '+file,sline
spawn,/sh,'wc -w '+file,sword
lines = long(sline)
words = long(sword)

n1 = words[0]/lines[0]
n2 = lines[0]

do_flat = keyword_set(flatten)
if do_flat then begin
    l = dindgen(n1>n2)
    ll = l*(l+1)/(2*!dpi)
endif

;help,n1,n2
if (n1 eq 1) then begin
    readcol, file, data, form='(d)', /silent
    data = reform(data, n2)
    if (do_flat) then data *= ll
endif else begin
    data = dblarr(n1, n2)
    openr, lun, file, /get_lun
    myline = dblarr(n1)
    for i=0, n2-1 do begin
        readf, lun, myline
        data[0, i] = myline
    endfor
    if (do_flat) then data *= ll#ll
endelse

return, data
end


