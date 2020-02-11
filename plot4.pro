close,/all
nx=10001

x1=dblarr(nx)
x2=dblarr(nx)
x3=dblarr(nx)
x4=dblarr(nx)

y1=dblarr(nx)
y2=dblarr(nx)
y3=dblarr(nx)
y4=dblarr(nx)

openr,15,'datos.dat'
for i=0,nx-1 do begin
 readf,15,a,b,c,d,e,f,g,h
x1(i)=a
y1(i)=b
x2(i)=c
y2(i)=d
x3(i)=e
y3(i)=f
x4(i)=g
y4(i)=h

endfor

 close, 15

set_plot, 'ps'					
loadct, 5
device, filename='4.ps', /color
loadct, 5

plot,x1,y1, back=255, col=0, thick=1, xtitle='', $
 ytitle='', chars=1, xr=[-2,2], yr=[-2,2]

oplot,x2,y2,col=15
oplot,x3,y3,col=30
oplot,x4,y4,col=45



device, /close					;cierro el device
set_plot, 'x'

end
