rgradrho = (rho(k)-rho(k-1))/dx
fgradrho = (rho(k+1)-rho(k))/dx
rgradv = (v(k)-v(k-1))/dx
fgradv = (v(k+1)-v(k))/dx
rgradp = (p(k)-p(k-1))/dx
fgradp = (p(k+1)-p(k))/dx
rgradeps = (eps(k)-eps(k-1))/dx
fgradeps = (eps(k+1)-eps(k))/dx



pro mhd04
; Computes and plots solutions to 1-D shock tube problem from Sod (1978)
; using finite difference methods, specifically the upwind (advection) scheme.
;
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;      Constants and initialization of variables    ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
npts=200			; commonly used...
dx = 1/200.d0 			; delta-x
total_time = 0.25d0
j = 0				; time step count
t = 0				; current time
dt = 0 				; time step 
; To satisfy Courant condition, the time step must be computed
; for each iteration, for now just initialize to zero
rhol=1.d0	; initial density left of the membrane
pl=1.d0		; initial pressure left of the membrane
rhor=0.125d0	; initial density right of the membrane
pr=0.1d0	; initial pressure right of the membrane
g=5/3.d0	; g(amma)
n=npts		; n columns, index k. each column represents a space grid point
m=3		; m rows, index i: i=0 for density, i=1 for dens*vel, i=2 for total energy density
;
; Define the functions to compute finite differences:
;
; "x" will hold actual lengths along the grid for final plotting.
x=dblarr(n) 
for k=0,n-1 do x(k)=(k+1)/200.	; Tube has unit length.
;
; The velocity is zero everywhere initially.
; The pressure is $p_r$ right of the membrane at x=1/2 and $p_l$ left of x=1/2. 
v = dblarr(n) 
p = dblarr(n)
for k= 0,(n-1)/2 do p(k) = pl
for k= (n-1)/2+1, n-1 do p(k) = pr
;
; "c" will hold values of \sqrt(v^2 + $c_{s}^2$) for computing a dt that satisfies Courant condition.
c=dblarr(n)
;
; row 0 (of a or b) is the density, "rho"
; row 1 is the momentum-density (velocity times density, "rho*v")
; row 2 is "epsilon" the total energy density= p/(g-1) + 0.5*rho*v*v
;
; intialize a.  for j (current time) = 0, v=0 so energy is just pressure over gamma-1
for k=0,(n-1)/2 do begin
  a(k,0) = rhol  
; a(k,1) = 0 by default
  a(k,2) = p(k) / (g-1)
endfor
for k=((n-1)/2+1),(n-1) do begin 	; Membrane at the halfway mark.
  a(k,0) = rhor  
; a(k,1) = 0 by default
  a(k,2) = p(k) / (g-1)
endfor
;
rho = a(*,0)  
rhov = a(*,1)
eps = a(*,2)
;
; Note epsilon is not to be confused with "u" the specific internal energy U/m from ideal gas law.
;
;
; Initialize conservation checks
M = 0
E = 0
; approximate integrals of density and energy density over length with discrete summation
for k=0,n-1 do begin
  M = M + dx*rho(k)
  E = E + dx*eps(k)
endfor
; store integrals for the tube at time zero, append later in the finite difference loop
msum = [M]
esum = [E]
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
;
; 
;
;
;
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;       Finite difference computation: Upwind Method     ;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
;
while (t le 0.25) do begin
  ;
  ;
  ;
  ; First we need to set the time step dt to satisfy the Courant condition
  ; so we need the sound speed squared at every grid point, which is the vector
  sqcs = g*p/rho
  ;
  ; Vector "c" holds dt values computed using the Courant condition
  ; for every grid point.
  for k=0,(n-1) do c(k) = dx / sqrt(v(k)*v(k) + sqcs(k))
  ;
  ; Now to find the minimum over the entire grid for this instant of time.
  ; Once we get the minimum, set dt equal to half the minimum just to be safe.
  min = c(0)
  for k=0,(n-1) do begin
    if (min gt c(k)) then min=c(k)
  endfor
  dt = 0.5*min
  ;
  ;
  ;
  ; Now we are ready to compute the evolution of the conditions in the tube
  ; forward through time.
  ;
  ; Finite differences:  the fact that the upwind piecewise function is 
  ; based on velocity sign is reflected in the conditional structure below.
  ;
  ;
  ; Continuity equation in advection form.  Density is advected.
  ; i=0 in index (k,i)
  for k=0,n-1 do begin
    ; Condition 1: no flow, everything stays the same
    ; Condition 2: v>0 and k=leftmost pt. No media behind it, so no flow from behind it.
    ; Condition 3: v<0 and k=rightmost pt. No media in front, so no flow through it from ahead.
    ; If none of these three are true, then k is v<0 or v>0, and furthermore the
    ; sign of v will pose no problem for the remaining combinations w/ endpts.
    if ( (v(k) eq 0) or ((v(k) gt 0) and (k eq 0)) or ((v(k) lt 0) and (k eq n-1)) ) $
    then b(k,0)=rho(k) else begin
      if (v(k) gt 0) then b(k,0) = rho(k) - v(k)*dt*(rho(k) - rho(k-1))/dx $
      else b(k,0) = rho(k) - v(k)*dt*(rho(k+1) - rho(k))/dx
    endelse
  endfor  
  ;
  ;
  ; Applying the given formulas for advection/upwind discretization naively, I found that
  ; the shock wouldn't propagate because after a time step, if v had changed, the 
  ; corresponding energy did not reflect that change, and it should since energy depends on v.
  ; So we have an interdependency problem, and we have to dig deeper into the advection differencing
  ; and pull out something more robust.
  ; If we consider the product rule applied to $\partial_t(\rho*v)$ and $\partial_x(\rho*v)$ then we 
  ; get new discretizations where density and velocity are "untangled" for lack of a better term.
  ; If we substitute those discretizations into the advection differencing for $\rho*v$ then 
  ; we get something very promising,
  ;
  ; Euler equation in advection form, density*velocity is advected.
  ; (Pressure is a "psuedo-advection" of density * sound speed.)
  ; i=1 in index (k,i)
  for k=0,n-1 do begin
    ;
    ; When the v = 0, pressure gradients can still exist, and they better
    ; exist or the shock won't propagate and there will never be any advection.
    ; In the Euler equation, when matter velocity is zero, we also have to take into 
    ; account the sound speed.  Sound speed can act analogously to advection via pressure 
    ; since $p = (1/\gamma) (\rho * c_s) * c_s$.
    ; (However, we don't take sound speed into account for the continuity equation in
    ; advection form because pressure is not in the continuity equation.)
    ; The choice of the direction of the pressure difference depends on the location.
    ; Left endpoint must have a forward difference, right end must have a reverse diff.
    ; Anywhere in between, it doesn't matter, but try an average of forward and reverse.
    ;
    ; v > 0, no problem w/ right endpoint, but for left endpt, no flow since nothing behind.
    ; v < 0, no problem w/ left endpoint, but for right endp, no flow since nothing in front.
    ;
    if (v(k) eq 0) then begin
      if (k eq 0) then b(k,1) = rhov(k) - dt*(p(k+1) - p(k))/dx
      if (k eq n-1) then b(k,1) = rhov(k) - dt*(p(k) - p(k-1))/dx
      if ((k gt 0) and (k lt n-1)) $
      then b(k,1) = rhov(k) - v(k)*dt*rhov(k)/dx - dt*(p(k+1) - p(k-1))/2/dx
    endif else begin
      if (v(k) gt 0) then begin
        if (k eq 0) then b(k,1) = rhov(k) - v(k)*dt*rhov(k)/dx - dt*(p(k+1) - p(k-1))/2/dx $
        else b(k,1) = rhov(k) - v(k)*dt*(rhov(k) - rhov(k-1))/dx - dt*(p(k) - p(k-1))/dx
      endif else begin  
        ; v(k) < 0 is true
        if (k eq n-1) then b(k,1) = rhov(k) + v(k)*dt*rhov(k)/dx + dt*p(k)/dx $
        else b(k,1) = rhov(k) - v(k)*dt*(rhov(k+1) - rhov(k))/dx - dt*(p(k+1) - p(k))/dx
      endelse
    endelse
    ;
    ;
  endfor  
  ;
  ;
  ;
  ; Energy Equation, both total energy density and pressure are advected
  ; since both total energy density and pressure gradients are multiplied by v.
  ; Therefore the conditional structure is just like density .
  ; energy density, i=2 in index (k,i)
  for k=0,n-1 do begin
    if ( (v(k) eq 0) or ((v(k) gt 0) and (k eq 0)) or ((v(k) lt 0) and (k eq n-1)) ) $
    then b(k,2)=eps(k) else begin
      if (v(k) gt 0) then b(k,2) = eps(k) - v(k)*dt*(eps(k) - eps(k-1))/dx $
      - v(k)*dt*(p(k) - p(k-1))/dx $
      else b(k,2) = eps(k) - v(k)*dt*(eps(k+1) - eps(k))/dx - v(k)*dt*(p(k+1) - p(k))/dx
    endelse
  endfor  
  ;
  ;
  ;
  ; Conservation checks - append the array that will hold the sum
  Mnext = 0
  Enext = 0
  for k=0,n-1 do begin
    Mnext = Mnext + dx*b(k,0)
    Enext = Enext + dx*b(k,2)
  endfor
  ; append summation array 
  msum = [msum, Mnext]
  esum = [esum, Enext]
  sums = [[msum],[esum]]
  ;
  ;
  ; Now set a, v, rho, p, and s for the next run, and current time j by time step.
  ; When j hits 0.25 or more, the while loop will end and the following will set the
  ; final values which we will plot:
  a = b
  rho = a(*,0)  
  rhov = a(*,1)
  eps = a(*,2)
  v = rhov / rho 
  p = (g-1) * (eps - 0.5*rhov*v)
  ; s = sackur-tetrode plus (p2/p1 -1)^3, a very tricky computation!
  j = j + 1
  t = t + dt
  ;
  ; stop 	; what the f is going on?
  ;
  ;
  ;
endwhile
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
;
; stop
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;       Plotting results        ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
set_plot,'ps'
;
device,filename='mhd04no2vel.ps'
plot, x,v, xtitle='x', ytitle='velocity';, psym=3
;
device,filename='mhd04no2dens.ps'
plot, x,rho, xtitle='x', ytitle='density';, psym=3
;
device,filename='mhd04no2press.ps'
plot, x,p, xtitle='x', ytitle='pressure';, psym=3
;
; device,filename='mhd04no2entropy.ps'
; plot, x,s, xtitle='x', ytitle='entropy';, psym=3
;
device,/close
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
;
end
