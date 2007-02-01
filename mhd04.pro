function fd, dir, var, m
  ; We're solving the Sod problem, and need finite differences, "fd"
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
      ; forward finite difference
        if dir eq 1 then diff = var(m+1) - var(m) 
      ;
      ; reverse finite difference
        if dir eq -1 then diff = var(m) - var(m-1)
      ;
      endelse
    endelse
    ;
  return, diff
end




function visc, direct, dens, vel, i
  ;
  ; To neglect artificial visc., change the following variable to 1.
  ; Now I don't have to erase every single call to viscosity computation
  ; in the main program to neglect viscosity for testing purposes.
  no_thanks = 1
  if (no_thanks eq 1) then begin
    Q = 0
    return, Q
  endif
  ;
  ; von Neumann and Richtmeyr Artificial Viscosity (from our pdf lecture notes)
  ;
  ; When viscosity is needed, it always occurs with specified direction of dv.
  ; Then each call to viscosity should pass the appropriate directional argument.
  ;
  ; dens and vel should be vectors of those respective quantities
  ; i indexes grid point where we want to know the viscosity.
  ;
  ; There are points (endpoints and adjacents) where we may not be able to
  ; compute Q because the velocity difference may not be computable at those points.
  ; No worries, the call to fd will take care of that, setting the velocity 
  ; difference 0, and then Q will be zero.  End result, no damping at tube edges.
  ;
  ; IDL doesn't distinguish between q and Q so had to name carefully.
  ;
  qconst=1.5  ; usually 0.05 < q < 2, as said in lecture notes
  ;
  ; Regain the velocity difference that the call to viscosity occured with
  dvel = fd(direct,vel,i)
  ;
  ; Now compute based on definition of Q
  if dvel lt 0 then Q = qconst^2*dens(i)*dvel^2 else Q = 0
  ;
  ;
  return, Q
end



function dvisc, direx, rhou, vee, indx
  ;
  ; The no viscosity option again.
  no_thanks = 1
  if (no_thanks eq 1) then begin
    dQ = 0
    return, dQ
  endif
  ;
  ;
  ; Same logical procedure as the "visc" function but computes a difference
  ; to match the direction of the velocity gradient occuring with the call to "dvisc"
  ;
  mesh = 200
  ;
  qconst=1.5  
  ;
  ; Q is still q^2*dens(k)*dvel(k)^2, but when differencing Q in a certain 
  ; direction, Qk - Qk-1, Qk+1 - Qk-1, or Qk+1 - Qk, you need the corresponding 
  ; pair of vel(k +/- 1,0), each vel. differenced in the specified direction. 
  ; Also, certain direction differences can't be computed near or at endpoints.
  ;
  if direx eq -1 then begin	; The reverse finite difference
    if indx lt 2  then dQ = 0 else begin
      Q = visc(direx,rhou,vee,indx)
      Ql = visc(direx,rhou,vee,indx-1)
      dQ = 0.5*(Q - Ql)
    endelse 
  endif
  ;
   if direx eq 1 then begin	; The forward finite difference
    if indx gt mesh-3  then dQ = 0 else begin
      Qr = visc(direx,rhou,vee,indx+1)
      Q = visc(direx,rhou,vee,indx)
      dQ = 0.5*(Qr - Q)
    endelse 
  endif
  ;   
  ;
  return, dQ
end




function find_waves, density
  ; Only the density has a discontinuity across every region boundary.
  ; So density will make a good marker for every wave position.
  ;
  ; The user should pass rho from the solution set.
  ;
  ; This function returns an array with the grid indices of the wave
  ; positions in row 0, and x values of the wave positions in row 1, with
  ; 3 columns, columns 0, 1, 2, are for rf, cd, sh, respectively.
  ;
  ;
  num = 200
  ;
  for grid=0,(num-1)/2 do begin
    while (density(grid+1) eq density(grid)) do begin
      distrf = (grid+1)/(num+0.d0)
      gridrf = grid+1
    endwhile
  endfor
  ;
  for grid=(num-1)/2,num-1 do begin
    while (density(grid+1) eq density(grid)) do begin
      distcd = (grid+1)/(num+0.d0)
      gridcd = grid+1
    endwhile
  endfor
  ;
  for grid=gridcd,num-1 do begin
    while (density(grid+1) eq density(grid)) do begin
      distsh = (grid+1)/(num+0.d0)
      gridsh = grid+1
    endwhile
  endfor
  ;
  wave_positions = [ [gridrf, gridcd, gridsh], [distrf, distcd, distsh] ]
  ;
  return, wave_positions
end




pro means_mhd04code
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
dx = 1/(n+0.d0) 	; delta-x
t = 0			; current time
dt = 0 			; time step, to be computed later to meet Courant condition
rhol=1.d0		; initial density left of the membrane at x=1/2
pl=1.d0			; initial pressure left of the membrane
rhor=0.125d0		; initial density right of the membrane
pr=0.1d0		; initial pressure right of the membrane
g=5/3.d0		; g(amma)
;
; The following are needed to compute specific entropy for an ideal gas
; according to the Sackur-Tetrode Equation. 
; (This is the only way I can see to get entropy w/ the information provided.)
;;
 ;
 kb = 1.38e-23
 plh = 6.626e-34
 ga = 5/3.d0
 mh = 1.673e-23
 mu = 29*mh 		; an average of Z for He and Xe
 ;
 ; Yeah i know in class we used H_2 and Xe but gamma=5/3 is for
 ; monatomic gas and H_2 ain't monatomic.  So pick the next 
 ; lightest thing, He.  I guess I could have assumed it was atomic
 ; hydrogen gas like in space but then 1-D isn't appropriate.
;;
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;   Initialize physical variables on the grid   ;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;  Finite difference computation: Advection differencing by Upwind Method  ;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
    fgrho = fd(1,rho,k)/dx
    fgv = fd(1,v,k)/dx
    fgrhov = fd(1,rhov,k)/dx
    fgp = fd(1,p,k)/dx
    fgrhoe = fd(1,rhoe,k)/dx
    ;
    ;
    ;
    if v(k) ge 0 then begin	; use reverse differences, put v=0 case here too since v>0 more likely than v<0
      ;
      nextrho(k) = rho(k) - v(k)*dt*rgrho - rho(k)*dt*rgv
      nextrhov(k) = rhov(k) - v(k)*dt*rgrhov - rhov(k)*dt*rgv $
      - dt*(rgp + dvisc(-1,rho,v,k)/dx)
      nextrhoe(k) = rhoe(k) - v(k)*dt*rgrhoe $
      - v(k)*dt*(rgp + dvisc(-1,rho,v,k)/dx) $
      - (rhoe(k)+p(k) + visc(-1,rho,v,k))*dt*rgv
      ;
    endif else begin		; use forward differences, v<0 not likely though
        ;
        nextrho(k) = rho(k) - v(k)*dt*fgrho - rho(k)*dt*fgv
        nextrhov(k) = rhov(k) - v(k)*dt*fgrhov - rhov(k)*dt*fgv $
	 - dt*(fgp + dvisc(1,rho,v,k)/dx)
        nextrhoe(k) = rhoe(k) - v(k)*dt*fgrhoe $
	- v(k)*dt*(fgp + dvisc(1,rho,v,k)/dx) $
        - (rhoe(k)+p(k) + visc(1,rho,v,k))*dt*fgv
	;
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
  ; 
  ; update other quantities for next run
  sqcs = g*p/rho
  c = dx / sqrt(v^2 + sqcs)
  ;
  ;
endwhile
;
;
;  Now that the evolution has completed, create entropy vector for plotting.
;  Entropy is set in a piecewise fashion.  But you need to know wave positions.
;
;myWaves = find_waves(rho)
;xsh = myWaves(1,2)
;xcd = myWaves(1,1)
;l = xsh - xcd			; Length of the post-shock region where entropy has jumped up.
;ksh = myWaves(0,2)
;kcd = myWaves(0,1)
;
; The entropy is calculated using the Sackur-Tetrode equation, the only way I can
; find to compute entropy (as opposed to just the change in entropy) using
; macroscopic variables instead of the partition function.  The change in entropy
; was given in our notes, and the change in entropy seems to follow from the
; Sackur-Tetrode, so, it all seems right.
; s = dblarr(n)
; for k=0,kcd do s(k)=5		; ...Xmas, the time of year for desserts like fudge.
; for k=ksh,n-1 do s(k)=17	; ...Lots of leftover fudge.
; for k=kcd,ksh do begin	; I trust my fudge factors, do you?
  ;
;  s(k) = 7 + kb*1e22/rho(k)/l*( 2.5 + 1.5*alog(p(k)*l/(g-1)) + alog(l) $
;  + alog( mu^4/rho(k)/l/plh^3*(4*!pi/3/rho(k)/l)^1.5 ) )
  ;
;endfor
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
;
;
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;       Plotting results        ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Plot analytical solution from #1 over these finite difference results
;
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
  px(k) = -0.3*((k+1)/(n+0.d0) - 1.5)^5
  rhx(k) = -0.5*((k+1)/(n+0.d0) - 1.5)^3
  vx(k) = 2.8*((k+1)/(n+0.d0) - 0.175)
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
;
;
; Plotting Numerical solution with analytical
;
loadct,40			; Graphing analytical and numerical in diff colors
dummy = dblarr(n)        	; Plots nothing, so we can get axes in black
;
set_plot,'ps'
;
device,filename='mhd04no2vel.ps'
plot, x,dummy, xtitle='x', ytitle='velocity'
oplot, x,v, color=60
oplot, x,vx, color=230, thick=3
;
device,filename='mhd04no2dens.ps'
plot, x,dummy, xtitle='x', ytitle='density'
oplot, x,rho, color=60
oplot, x,rhx, color=230, thick=3
;
device,filename='mhd04no2press.ps'
plot, x,dummy, xtitle='x', ytitle='pressure'
oplot, x,p, color=60
oplot, x,px, color=230, thick=3
;
;device,filename='mhd04no2entropy.ps'
;plot, x,dummy, xtitle='x', ytitle='entropy'
;oplot, x,s, color=60
;oplot, x,sx, color=230, thick=3
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

    
