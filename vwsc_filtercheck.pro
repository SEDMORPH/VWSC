;;*** AIM: check filters sent from Chris against the ones I have been
;;using. Chris's Spitzer filters are much smoother, but other
;;than that there is very little difference. 
 
;; email from 6th March 2014
;; "In principle, u, B, V, i (called ip), z (zp), ch1 and ch2 should be the same as for the UDS (same telescope and instrument), or at least as close as makes no difference."

PRO VWSC_FILTERCHECK

;;-- my filters
name1 = [ 'Subaru/suprime_B.txt','Subaru/suprime_V.txt', 'Subaru/suprime_i.txt','Subaru/suprime_z.txt','VISTA/vircam_Y.txt','SPITZER/3p6_spitzer.txt','SPITZER/4p5_spitzer.txt']

;;-- from Chris Simpson
name2 = [ 'UDSliv/filt_B_vw.dat','UDSliv/filt_V_vw.dat','UDSliv/filt_ip_vw.dat','UDSliv/filter_z_vw.dat','UDSliv/filt_Y_vw.dat','UDSliv/filt_ch1_vw.dat','UDSliv/filt_ch2_vw.dat']

name3 = ['UDSliv/filt_B.dat','UDSliv/filt_V.dat','UDSliv/filt_ip.dat','UDSliv/filter_z.dat','UDSliv/filt_Y.dat','UDSliv/filt_ch1.dat','UDSliv/filt_ch2.dat']

ps1,'~/home/idl/pro/VWSC/filtercheck.ps'
for i=0,n_elements(name1)-1 do begin
   readcol,!dataDIR+'FILTERS/'+name1[i], lam1,fil1,form='(F,F)'
   readcol, !dataDIR+'FILTERS/'+name2[i], lam2,fil2,form='(F,F)'
   readcol, !dataDIR+'FILTERS/'+name3[i], lam3,fil3,form='(X,F,F)'

   plot, lam1,fil1/max(fil1),/xs,title=name1[i] + ' ' + name2[i]      
   oplot, lam3,fil3/max(fil3),color=cgcolor('red')
   oplot, lam2,fil2/max(fil2),color=cgcolor('blue')

endfor
ps2

END
