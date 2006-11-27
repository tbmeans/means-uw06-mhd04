function fd, dir, var, m
  ; We're solving the Sod problem, and need finite differences
  ; Takes a vector 'var' holding values of a variable over the grid
  ; The dir parameter should be 0,1, or -1, to determine difference direction
  ; m should be the index of grid pt around which a difference is needed
  ; and where such a difference is possible (sometimes not possible at endpts).
  ; This function will return differences of zero when:
  ;  - it is asked to give anything but a forward difference at left endpoint.
  ;  - it is asked to give a anything but a reverse difference at right endpt.
  ;  The ultimate result of these exceptions will be that the variable doesn't
  ;  change for the current time step.
    ;
    npts = 200  ; This is set for the main program later as well
    ;
    if ((m eq 0) and (dir ne 1)) then diff = 0 else begin
      if ((m eq npts-1) and (dir ne -1)) then diff = 0 else begin
      ;
      ; standard finite difference
        if dir eq 0 then diff = (var(m+1)-var(m-1))/2
      ;
      ; forward finite difference
        if dir eq 1 then diff = var(m+1) - var(m) 
      ;
      ; reverse finite difference
        if dir eq -1 then diff = var(m) - var(m-1)
      ;
      endelse
    endelse
    ;
  return diff
end




function visc, dens, vel, i
  ; von Neumann and Richtmeyr Artificial Viscosity (from our pdf lecture notes)
  ;
  ; dens and vel should be vectors of those respective quantities
  ; i indexes grid point where we want to know the viscosity
  ; The user should never pass i=0,1,npts-2, nor npts-1
  ; npts should be the number of gridpoints, xstep the grid step.
  ;
  ;
  qconst=1.5  ; usually 0.05 < q < 2, as said in lecture notes
  ; I need a lot of damping, based on my solution w/o artificial visc.
  ; IDL doesn't distinguish between q and Q so had to name carefully.
  ;
  ; By definition of Q, we can't calculate Q at the endpts nor at the point
  ; adjacent to each endpoint, so we'll say no damping at tube edges.
  ; Won't need it there anyway.
  ; Then decide on forward or reverse velocity gradients based on sign of v
  ; That will make this consistent with upwind scheme.
  ;
  dvfwd = fd(1,vel,i)		
  dvrev = fd(-1,vel,i)			; v diffs select piecewise Q
  dv = fd(0,vel,i)
  Q = qconst^2*dens(i)*dv^2			; Q at given point i
  Qr = qconst^2*dens(i+1)*(fd(0,vel,i+1))^2	; Q at i+1
  Ql = qconst^2*dens(i-1)*(fd(0,vel,i-1))^2	; Q at i-1
  ;
  if vel(i) lt 0 then begin
    if dvfwd lt 0 then dQ = Qr - Q else begin
      Q = 0
      dQ = 0
      ; Well... that's what the pdf lecture notes said...
    endelse
  endif else begin	
    if vel(i) gt 0 then begin
      if dvrev lt 0 then dQ = Q - Ql else begin		
	Q = 0
        dQ = 0
      endelse
    endif else begin 	; only v=0 is left
      if dv lt 0 then dQ = 0.5*(Qr - Ql) else begin
        Q = 0 
	dQ = 0		
      endelse
    endelse
  endelse
  ;
  ;
  ; The kill switch.  To neglect artificial viscosity,
  ; change the following variable to 1.  Now I don't have to
  ; erase every single call to viscosity value or computation
  ; in the main program to neglect viscosity for testing purposes.
  ;
  no_thanks = 0
  ;
  if (no_thanks eq 1) then begin
    Q = 0
    dQ = 0
  endif
  ;
  ;
  Qset = [Q, dQ]
  return Qset
end
  



function st, ps, rh, xes
  ; computing sackur-tetrode specific entropy
  ; soln should be the solution matrix with rho, rho*v, rho*e
  ; the user should pass a set of x and k for wave positions.
  ;
  kb = 1.38e-23
  plh = 6.626e-34
  ga = 5/3.d0
  mh = 1.673e-23
  mu = 29*mh 		; an average of Z for He and Xe
  ;
  ; Yeah i know in class we used H_2 and Xe but gamma=5/3 is for
  ; monatomic gas and H_2 ain't monatomic.  So pick the next 
  ; lightest thing, He.
  ;
  nizzle=200
  ste = dblarr(nizzle)
  for a=0,kcd do ste(a)=5
  for a=ksh,nizzle-1 do ste(a)=17
  l = xsh - xcd
  for b=kcd,ksh do begin	; just trust the fudge factors...
    ste(b) = 7 + k*1e22/rh(b)/l*( 2.5 + 1.5*alog(ps(b)*l/(g-1)) + alog(l) $
    + alog( mu^4/rh(b)/l/plh^3*(4*!pi/3/rho(b)/l)^1.5 ) )
  endfor
  ;
  return ste
end




function findxs, density
  ; Only the density has a discontinuity across every region boundary
  ; So density will make a good marker for every wave position
  ; The user should pass rho from the solution set.
  ;
  num = 200
  ;
  for k=0,(num-1)/2 do begin
    while (density(k+1) eq density(k)) do begin
      xrf = (k+1)/(num+0.d0)
      krf = k+1
    endwhile
  endfor
  ;
  for k=(num-1)/2,num-1 do begin
    while (density(k+1) eq density(k)) do begin
      xcd = (k+1)/(num+0.d0)
      kcd = k+1
    endwhile
  endfor
  ;
  for k=kcd,num-1 do begin
    while (density(k+1) eq density(k)) do begin
      xsh = (k+1)/(num+0.d0)
      ksh = k+1
    endwhile
  endfor
  ;
  xkset = [[krf, kcd, ksh],[xrf, xcd, xsh]]
  return xkset
end




pro mhd04
; Computes and plots solutions to 1-D shock tube problem from Sod (1978)
; using finite difference methods, specifically the upwind (advection) scheme.
;
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;      Constants     ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
n=200			; commonly used number of grid points
dx = 1/200.d0 		; delta-x
t = 0			; current time
dt = 0 			; time step, to be computed later to meet Courant condition
rhol=1.d0		; initial density left of the membrane at x=1/2
pl=1.d0			; initial pressure left of the membrane
rhor=0.125d0		; initial density right of the membrane
pr=0.1d0		; initial pressure right of the membrane
g=5/3.d0		; g(amma)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;   The simulation grid   ;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
x=dblarr(n)     		; holds distances x along the grid for final plotting.
for k=0,n-1 do x(k)=(k+1)/(n+0.d0)	; Tube has unit length.
;
v = dblarr(n) 		; Velocity is zero everywhere initially.
p = dblarr(n)		; Pressure is $p_r$ right of x=1/2 and $p_l$ left of x=1/2. 
rho = dblarr(n)		; Density is $\rho_r$ right of x=1/2 and $\rho_l$ left of x=1/2.
;
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;   Intialize vel., density and pressure.   ;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
for k= 0,(n-1)/2 do begin	; Grid pts 0 to 99 for x 0 to 1/2
  p(k) = pl
  rho(k) = rhol
endfor
for k= (n-1)/2+1, n-1 do begin	; Grid pts 100 to 199 for x 1/2 to 1
  p(k) = pr
  rho(k) = rhor
endfor  
;
sqcs = g*p/rho		; The sound speed squared
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;  Initialize other physical variables  ;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
rhov = rho*v 			; Momentum-density.  
rhoe = p/(g-1) + 0.5*rhov*v	; density*total specific E = total E density
c = dx / sqrt(v^2 + sqcs)	
;
; Hereafter, we treat rho*v as an inseparable product, a single variable.
; "c" is for values of $\sqrt(v^2 + c_{s}^2)$ for dt to meet Courant condition.  
;
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;   Initialize conservation checks   ;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
M = 0
E = 0
; approximate integrals of density and energy density over length with discrete summation
for k=0,n-1 do begin
  M = M + dx*rho(k)
  E = E + dx*rhoe(k)
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
  ; where dt < min( dx / sqrt(v*v + cs*cs) )
  ;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;   Minimization routine   ;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  min = c(0)
  for k=0,(n-1) do begin
    if (min gt c(k)) then min=c(k)
  endfor
  dt = 0.01*min
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;  
  ;
  ;  
  ; Values after this timestep
  nextrho = dblarr(n)
  nextrhov = dblarr(n)
  nextrhoe = dblarr(n) 
  ;
  for k=0,n-1 do begin  ; Evolving the system in time for every grid point.
    ;
    ; Prepare discrete gradients
    ;
    rgrho = fd(-1,rho,k)/dx
    rgv = fd(-1,v,k)/dx
    rgrhov = fd(-1,rhov,k)/dx
    rgp = fd(-1,p,k)/dx
    rgrhoe = fd(-1,rhoe,k)/dx
    ;
     grho = fd(0,rho,k)/dx
     gv = fd(0,v,k)/dx
     grhov = fd(0,rhov,k)/dx
     gp = fd(0,p,k)/dx
     grhoe = fd(0,rhoe,k)/dx
    ;
    fgrho = fd(1,rho,k)/dx
    fgv = fd(1,v,k)/dx
    fgrhov = fd(1,rhov,k)/dx
    fgp = fd(1,p,k)/dx
    fgrhoe = fd(1,rhoe,k)/dx
    ;
    ;
    ; Prepare artificial viscosity.  fd and visc functions ensure proper
    selection of viscosity difference direction.
    ;
    vsc_set = visc(rho,v,k)
    vsc = vsc_set(0)
    gvsc = vsc_set(1) / dx
    ;
    ; fd and visc functions ensure proper selection of viscosity 
    ; difference direction.
    ;
    ;
    if v(k) gt 0 then begin	; use reverse differences
      ;
      nextrho(k) = rho(k) - v(k)*dt*rgrho - rho(k)*dt*rgv
      nextrhov(k) = rhov(k) - v(k)*dt*rgrhov - rhov(k)*dt*rgv - dt*(rgp+gvsc)
      nextrhoe(k) = rhoe(k) - v(k)*dt*rgrhoe - v(k)*dt*(rgp+gvsc) $
      - (rhoe(k)+p(k)+vsc)*dt*rgv
      ;
    endif else begin
      ;
      if v(k) lt 0 then begin   ; use forward gradients
        nextrho(k) = rho(k) - v(k)*dt*fgrho - rho(k)*dt*fgv
        nextrhov(k) = rhov(k) - v(k)*dt*fgrhov - rhov(k)*dt*fgv - dt*(fgp+gvsc)
        nextrhoe(k) = rhoe(k) - v(k)*dt*fgrhoe - v(k)*dt*(fgp+gvsc) $
        - (rhoe(k)+p(k)+vsc)*dt*fgv
      ;
      endif else begin	; v=0: average the forward and reverse gradients
        ;
	nextrho(k) = rho(k) - v(k)*dt*grho - rho(k)*dt*gv
        nextrhov(k) = rhov(k) - v(k)*dt*grhov - rhov(k)*dt*gv $
	- dt*(gp+gsvc)
        nextrhoe(k) = rhoe(k) - v(k)*dt*grhoe - v(k)*dt*(gp+gvsc) $
	- (rhoe(k)+p(k)+vsc)*dt*gv
      endelse
    endelse
    ;
    ;
  endfor
  ;
  ;
  ;
  ; Conservation checks - append the array that will hold the sum
  Mnext = 0
  Enext = 0
  for k=0,n-1 do begin
    Mnext = Mnext + dx*nextrho(k)
    Enext = Enext + dx*nextrhoe(k)
  endfor
  ; append summation array 
  msum = [msum, Mnext]
  esum = [esum, Enext]
  sums = [[msum],[esum]]
  ;
  ;
  ; Now set quantities for the next run, and current time t by time step dt.
  ; When t hits 0.25 or more, the while loop will end and the following will set the
  ; final values which we will plot:
  t = t + dt
  rho = nextrho
  rhov = nextrhov
  rhoe = nextrhoe
  v = nextrhov / nextrho 
  p = (g-1) * (rhoe - 0.5*rhov*v)
  wavex = findxs(rho)
  s = st(p, rho, wavex)  ; entropy defined piecewise based on wave positions
  ;
  ; update other quantities for next run
  sqcs = g*p/rho
  c = dx / sqrt(v^2 + sqcs)
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
; Plot analytical solution from #1 over these finite difference results
;
; Analytical solution
;
;  prf, rhorf, vrf, I don't know how to get the curvy things
;
akrf = 34
akcd = 144
aksh = 174
px = dblarr(n)
rhx = dblarr(n)
vx = dblarr(n)
sx = dblarr(n)
;

for k=0,akrf do begin
  px(k) = pl
  rhx(k) = rhol
  vx(k) = 0
  sx(k) = 5
endfor

for k=akrf+1,100 do begin
  px(k) = 0 ; some crazy shit i haven't got yet
  rhx(k) = 0 ; some crazy shit i haven't got yet 
  vx(k) = 0 ; supposed to be a diagnonal line
  sx(k) = 5
endfor

for k=101,akcd do begin
  px(k) = 0.3
  rhx(k) = 0.5
  vx(k) = 0.9
  sx(k) = 5
endfor

for k=akcd+1,aksh do begin
  px(k) = 0.3
  rhx(k) = 0.2
  vx(k) = 0.9
  sx(k) = 13
endfor

for k=akcd+1,aksh do begin
  px(k) = 0.3
  rhx(k) = 0.2
  vx(k) = 0.9
  sx(k) = 13
endfor

for k=aksh+1,n-1 do begin
  px(k) = 0.1
  rhx(k) = 0.125
  vx(k) = 0
  sx(k) = 17
endfor

;
set_plot,'ps'
;
device,filename='mhd04no2vel.ps'
plot, x,v, xtitle='x', ytitle='velocity', psym=4
oplot , psym=0
;
device,filename='mhd04no2dens.ps'
plot, x,rho, xtitle='x', ytitle='density'
;
device,filename='mhd04no2press.ps'
plot, x,p, xtitle='x', ytitle='pressure'
;
device,filename='mhd04no2entropy.ps'
plot, x,s, xtitle='x', ytitle='entropy'
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

    
