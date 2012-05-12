function solarv, timespec, x, y, hg=hg, rotmodel=rotmodel, $
                 status=status, verbose=verbose

  cs = 'xy'
  if keyword_set (hg) then cs = 'lola'
  if not keyword_set (rotmodel) then rotmodel = 'fixed'
  
  cmd = 'solarv -m ' + rotmodel + ' ' + str(timespec) $
        + ' ' + str(cs) + ' ' + str(x) + ' ' + str(y)

  if keyword_set (verbose) then $
     print, "calling cmd='" + str(cmd) + "'"

  spawn, cmd, result, exit_status=status

  if keyword_set (verbose) then $
     print, "result=", result

  if status eq 0 then begin
     fields = strsplit (result[1], ' ', /extract)
  endif else begin
     fields = dblarr (8)
  endelse
     
  res = { jd:  fields[0], $
          mjd: fields[0] - 2400000.5d, $
          b0:  fields[1], $
          l0:  fields[3], $
          vlos: fields[4], $
          dist: fields[5], $
          vlos_c: fields[6], $
          dist_c: fields[7] $
        }

  return, res
end
