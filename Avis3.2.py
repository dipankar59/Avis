# Avis3.2.py
# Syntax update for compatibility with both python 2.7 and python 3
#
# Mandar Hulsurkar and Dipankar Bhattacharya 19 May 2021
#
# Avis3.1.py
#
# Code to compute Astrosat orbital position, velocity and source 
# visibility as a function of time.
#
# Requires a Two-Line Element (TLE) file for an epoch close to
# the time at which the results are desired.
#
# Uses SGP4 propagation
#
# Incorporates a tunable trapezoidal SAA exclusion zone
#
# Outputs text files containing the table of parameters as a 
# function of time, for each source in the file "sourcelist.txt"
#
# Settable parameters supplied through an input file Avis3_config.txt
#
# Output visibility file contains instrument-wise total visibility
# for each orbit.
#
# A more detailed orbitview file can be printed if desired, which
# contains the state vectors and other parameters at every sample
# time (typically 1 sec to 1 min, settable through the config file).
#
# Dipankar Bhattacharya June 2014; last update Sep 2016
# Update from Avis3.py: negative roll constraint applied to all 
#                       instruments.  March 2017 
# 
#--------------------------------------------------------------

import os
from math import *
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4

print ("")
print ("## Astrosat Visibility Calculator ##")

# Default values
# Start Date, time  and Span over which orbit is to be computed
date=1,10,2016
ut0=0,0,0.0
SpanDays=10.0
SecInterval=1.0

# SAA model
SAAlatRef=-6.0
SAAlongMin=-110.0
SAAlongMax=0.0
SAAslopeLeft=3.5
SAAslopeRight=0.0

# View constraints
SunAvoid=65.0
MoonAvoid=15.0
LimbAvoid=12.0
RamAvoid=12.0
NegRollSunAvoid=30.0

# Output control
PrintOrb=1
PrintVis=1

# Read settable parameters from the config file Avis3_config.txt
cfgname="Avis3_config.txt"
if os.path.isfile(cfgname):
  configfile=open(cfgname,'r')
  for line in configfile:
    line.strip()
    columns=line.split()
    if (columns[0]=="StartDate"):
       date=int(columns[1]),int(columns[2]),int(columns[3])

    if (columns[0]=="StartUT"):
       ut0=int(columns[1]),int(columns[2]),float(columns[3])

    if (columns[0]=="SpanDays"):
       SpanDays=float(columns[1])

    if (columns[0]=="SecInterval"):
       SecInterval=float(columns[1])

    if (columns[0]=="SAAlatRef"):
       SAAlatRef=float(columns[1])

    if (columns[0]=="SAAlongMin"):
       SAAlongMin=float(columns[1])

    if (columns[0]=="SAAlongMax"):
       SAAlongMax=float(columns[1])

    if (columns[0]=="SAAslopeLeft"):
       SAAslopeLeft=float(columns[1])

    if (columns[0]=="SAAslopeRight"):
       SAAslopeRight=float(columns[1])

    if (columns[0]=="SunAvoid"):
       SunAvoid=float(columns[1])

    if (columns[0]=="MoonAvoid"):
       MoonAvoid=float(columns[1])

    if (columns[0]=="RamAvoid"):
       RamAvoid=float(columns[1])

    if (columns[0]=="LimbAvoid"):
       LimbAvoid=float(columns[1])

    if (columns[0]=="NegRollSunAvoid"):
       NegRollSunAvoid=float(columns[1])

    if (columns[0]=="PrintOrb"):
       PrintOrb=float(columns[1])

    if (columns[0]=="PrintVis"):
       PrintVis=float(columns[1])

  configfile.close()
else:
  print ("config file not found. Using default")

SpanSecs=SpanDays*86400.0

# Satellite orbital elements
tleline=['\n','\n', '\n']
tlename="Astrosat.tle"
if os.path.isfile(tlename):
  tlefile=open(tlename,'r')
  i=0
  for line in tlefile:
    line.strip()
    tleline[i]=line
    i+=1
  tlefile.close()
else:
  print ("TLE file not found. Using default")
  tleline[0]='Astrosat'
  tleline[1]=('1 40930U 15052A   16254.67150853  .00001010  00000-0  53365-4 0  9990')
  tleline[2]=('2 40930   5.9955  66.8928 0008517 253.8813 106.0572 14.76026389 51587')

dr=pi/180.
Rearth=6378.135
satellite=twoline2rv(tleline[1], tleline[2], wgs72)
Porbsat=120.0*pi/satellite.no
epochsat=satellite.jdsatepoch-2451543.5

def refday(date):
    dd,mm,yr=date
    obsday=367.0*yr-(7*(yr+((mm+9)//12)))//4+(275*mm)//9+dd-730530
    return obsday

def sunmoon(obsday):

# Sun elements
    ws=282.9404+(4.70935e-5)*obsday
    es=0.016709-1.151e-9*obsday
    Ms=356.0470+0.9856002585*obsday
    ob=23.4393-3.563e-7*obsday

# Moon elements
    N=125.1228-0.0529538083*obsday
    incl=5.1454
    wm=318.0634+0.1643573223*obsday
    am=60.2666
    em=0.054900
    Mm=115.3654+13.0649929509*obsday

    Ls=ws+Ms
    E0s=Ms+es*sin(dr*Ms)*(1.0+es*cos(dr*Ms))/dr

    x=cos(dr*E0s)-es
    y=sin(dr*E0s)*sqrt(1.0-es*es)
    r=sqrt(x*x+y*y)
    v=(atan2(y,x)/dr)%360

    lon=v+ws
    x=r*cos(dr*lon)
    y=r*sin(dr*lon)

    xsun=x
    ysun=y*cos(dr*ob)
    zsun=y*sin(dr*ob)
    rsun=sqrt(xsun*xsun+ysun*ysun+zsun*zsun)

    E0m=Mm+em*sin(Mm*dr)*(1.0+em*cos(Mm*dr))/dr
    Ediff=1.0
    while Ediff > 0.005:
        E1m=E0m-(E0m-em*sin(E0m*dr)/dr-Mm)/(1.0-em*cos(E0m*dr))
        Ediff=abs(E0m-E1m)
        E0m=E1m

    x=am*(cos(E0m*dr)-em)
    y=sin(E0m*dr)*am*sqrt(1.0-em*em)
    r=sqrt(x*x+y*y)
    v=(atan2(y,x)/dr)%360

    Na=N*dr
    vwa=(v+wm)*dr
    ia=incl*dr
    xec=r*(cos(Na)*cos(vwa)-sin(Na)*sin(vwa)*cos(ia))
    yec=r*(sin(Na)*cos(vwa)+cos(Na)*sin(vwa)*cos(ia))
    zec=r*sin(vwa)*sin(ia)
    lon=(atan2(yec,xec)/dr)%360
    lat=asin(zec/r)/dr
    
# Moon perturbations
    Lm=N+wm+Mm
    D=Lm-Ls
    F=Lm-N
    
    dlon=-1.274*sin(dr*(Mm-2.*D))+0.658*sin(2.*D*dr)-0.186*sin(Ms*dr)
    dlon=dlon-0.059*sin(2.*dr*(Mm-D))-0.057*sin(dr*(Mm-2*D+Ms))
    dlon=dlon+0.053*sin(dr*(Mm+2*D))+0.046*sin(dr*(2*D-Ms))
    dlon=dlon+0.041*sin(dr*(Mm-Ms))-0.035*sin(dr*D)-0.031*sin(dr*(Mm+Ms))
    dlon=dlon-0.015*sin(2*dr*(F-D))+0.011*sin(dr*(Mm-4.*D))

    dlat=-0.173*sin(dr*(F-2*D))-0.055*sin(dr*(Mm-F-2*D))
    dlat=dlat-0.046*sin(dr*(Mm+F-2*D))+0.033*sin(dr*(F+2*D))
    dlat=dlat+0.017*sin(dr*(2*Mm+F))
    
    delr=-0.58*cos(dr*(Mm-2*D))-0.46*cos(2*dr*D)

    lon=lon+dlon
    lat=lat+dlat
    r=r+delr
    
    rmoon=r*Rearth
    zec=rmoon*sin(dr*lat)
    yec=rmoon*cos(dr*lat)*sin(dr*lon)
    xec=rmoon*cos(dr*lat)*cos(dr*lon)

    xmoon=xec
    ymoon=yec*cos(ob*dr)-zec*sin(ob*dr)
    zmoon=yec*sin(ob*dr)+zec*cos(ob*dr)

    smpos=xsun,ysun,zsun,xmoon,ymoon,zmoon
    return(smpos)

def satpos(obsday,phi):
    dmphi=(phi*dr)/satellite.no
    minutes=(obsday-epochsat)*1440.0+dmphi
    sp,sv=sgp4(satellite,minutes)
    xs=sp[0]
    ys=sp[1]
    zs=sp[2]
    vx=sv[0]
    vy=sv[1]
    vz=sv[2]
    pvsat=xs,ys,zs,vx,vy,vz
    return(pvsat)

def dmsconv(deg,min,sec):
  if deg < 0:
     sign=-1
  else:
     sign=+1
  absdeg=sign*deg
  degf=absdeg+(float(min)+sec/60.0)/60.0
  degf=sign*degf
  return degf

def ParseInputLine(line):
  indata=line.split()
  name=indata[0]
  rastr=indata[1]
  decstr=indata[2]
  delim=":"
  if delim in rastr:
     rra=rastr.split(delim)
     hour=int(rra[0])
     min=int(rra[1])
     if len(rra) < 3 :
        sec=0.0
     else:
        sec=float(rra[2])
     RA2000=15.0*dmsconv(hour,min,sec)
  else:
     RA2000=float(rastr)
  if delim in decstr:
     rdec=decstr.split(delim)
     deg=int(rdec[0])
     min=int(rdec[1])
     if len(rdec) < 3 :
        sec=0.0
     else:
        sec=float(rdec[2])
     Dec2000=dmsconv(deg,min,sec)
  else:
     Dec2000=float(decstr)
  return(name,RA2000,Dec2000)

dayref=refday(date)
hh,mm,ss=ut0
dayref=dayref+float(hh)/24.0+float(mm)/1440.0+ss/86400.0
mjdref=dayref+51543.0

infile=open('sourcelist.txt','r')

for line in infile:
  line.strip()
  target,srcRA2000,srcDec2000=ParseInputLine(line)

# Precession coefficients
  fRAdeg=1.14077e-5*(3.075+1.336*sin(srcRA2000*dr)*tan(srcDec2000*dr))
  fDecdeg=1.52407e-5*cos(srcRA2000*dr)

# Prepare output files
  outname=target+"_orb.txt"
  visname=target+"_avis.txt"
  print ("%15s %9.4f %9.4f : " % (target,srcRA2000,srcDec2000))
  if (PrintOrb == 1):
     #outfile=open(outname,'wb')
     outfile=open(outname,'w')
     print ("   writing orbitview  file %s" % (outname))

  if (PrintVis == 1):
     #visfile=open(visname,'wb')
     visfile=open(visname,'w')
     print ("   writing visibility file %s" % (visname))

  if (PrintOrb == 1):
     outfile.write('# Astrosat Mission. Orbital Period %8.2f s\n' %(Porbsat))
     outfile.write("# Orbit and Elevation for %s " %(target))
     outfile.write("RA2000: %8.4f deg Dec2000: %8.4f deg\n" %(srcRA2000,srcDec2000))
     outfile.write("# Starting on "+str(date[0])+"/"+str(date[1])+"/"+str(date[2]))  
     outfile.write(" "+str(ut0[0])+":"+str(ut0[1])+":"+str(ut0[2])+"UT ")
     outfile.write(": MJD  %13.7f\n" %(mjdref))
     outfile.write("# Constraints: Sun %5.1f Moon %4.1f Ram %4.1f Limb %4.1f NegRollSun %5.1f\n" %(SunAvoid,MoonAvoid,RamAvoid,LimbAvoid,NegRollSunAvoid))
     outfile.write("# SAA Trapz LatRef %5.1f Long Min %6.1f Max %6.1f Slope L %5.2f R %5.2f\n" %(SAAlatRef,SAAlongMin,SAAlongMax,SAAslopeLeft,SAAslopeRight))
     outfile.write("# SAAflag: 1=in, 0=out   Nightflag: 1=night, 0=day   Visflag: 1=source observable, 0=not observable\n")
     outfile.write("# System UT,X,Y,Z,Vx,Vy,Vz,longitude,latitude,SunAngle,MoonAngle,Ramangle,Elevation,SAAflag,Nightflag,Visflag\n")
     outfile.write("#-----------------------------------------------------------------------------------------------------------------------------\n")

  if (PrintVis == 1):
     visfile.write('# Astrosat Mission. Orbital Period %8.2f s\n' %(Porbsat))
     visfile.write("# Visibility for %s " %(target))
     visfile.write("RA2000: %8.4f deg Dec2000: %8.4f deg\n" %(srcRA2000,srcDec2000))
     visfile.write("# Starting on "+str(date[0])+"/"+str(date[1])+"/"+str(date[2]))  
     visfile.write(" "+str(ut0[0])+":"+str(ut0[1])+":"+str(ut0[2])+"UT ")
     visfile.write(": MJD  %13.7f\n" %(mjdref))
     visfile.write("# Constraints: Sun %5.1f Moon %4.1f Ram %4.1f Limb %4.1f NegRollSun %5.1f\n" %(SunAvoid,MoonAvoid,RamAvoid,LimbAvoid,NegRollSunAvoid))
     visfile.write("# SAA Trapz LatRef %5.1f Long Min %6.1f Max %6.1f Slope L %5.2f R %5.2f\n" %(SAAlatRef,SAAlongMin,SAAlongMax,SAAslopeLeft,SAAslopeRight))
     visfile.write("# St.MET, elapsed days, MJD, min Ram, Sun ang, Moon ang, vis (s): LXP/CZT, SXT, UVT\n")
     visfile.write("#----------------------------------------------------------------------------------\n")

# Compute state vectors at each sample time, orbit-wise visibility and output results
  orbnum=0
  newSAA=0
  newNight=0
  SAAbegin=0.0
  SAAend=0.0
  Nightbegin=0.0
  Nightend=0.0
  dphase=float(SecInterval)/Porbsat
  orbphase=1.0-dphase
  elapsed=-SecInterval
  while elapsed < SpanSecs:
    elapsed=elapsed+SecInterval
    elapsed_days=float(elapsed)/86400.0
    obsday=dayref+elapsed_days
    mjd=obsday+51543.0
    orbphase=orbphase+dphase
    if orbphase >= 1.0:
       if orbnum > 0:
          if minram < RamAvoid:
             lxczvis=0.0
             sxtvis=0.0
             uvtvis=0.0
          uvtvis=uvtvis-370.0
          if ((SAAbegin>Nightbegin) and (SAAend<Nightend)):
             uvtvis=uvtvis-220.0
          if (uvtvis < 150.0):
             uvtvis=0.0
          if (PrintVis == 1):
             visfile.write('%12.2f %12.7f %13.7f %5.1f %5.1f %5.1f   %7.1f %7.1f %7.1f\n' % (orbT0,elapsed_days,mjd,minram,sunangle,moonangle,lxczvis,sxtvis,uvtvis))
       orbphase=orbphase-1.0
       orbT0=(obsday-3654.0)*86400.0-orbphase*Porbsat
       lxczvis=0.0
       sxtvis=0.0
       uvtvis=0.0
       minram=90.0
       orbnum=orbnum+1

# Precess source coordinates to those of date
    srcRA=srcRA2000+fRAdeg*(obsday-1)
    srcDec=srcDec2000+fDecdeg*(obsday-1)

# Greenwich Mean Sidereal Time
    GMST=(280.46061837+360.98564736629*(obsday-1))%360

    xsrc=cos(srcRA*dr)*cos(srcDec*dr)
    ysrc=sin(srcRA*dr)*cos(srcDec*dr)
    zsrc=sin(srcDec*dr)

    xsun,ysun,zsun,xmoon,ymoon,zmoon=sunmoon(obsday)
    rsun=sqrt(xsun*xsun+ysun*ysun+zsun*zsun)

    sunangle=(acos((xsrc*xsun+ysrc*ysun+zsrc*zsun)/rsun)/dr)%360

    phi=0.0
    xsat,ysat,zsat,vx,vy,vz=satpos(obsday,phi)
    satdist=sqrt(xsat*xsat+ysat*ysat+zsat*zsat)
    satspeed=sqrt(vx*vx+vy*vy+vz*vz)

    satlat=(asin(zsat/satdist))/dr
    satlong=((atan2(ysat,xsat))/dr)%360-GMST
    if satlong < -180.0 :
       satlong=satlong+360.0
    if satlong > 180.0 :
       satlong=satlong-360.0
    
    SAAshift1=SAAslopeLeft*(satlat-SAAlatRef)
    SAAshift2=SAAslopeRight*(satlat-SAAlatRef)

    SAAflag=0
    if ((satlong > (SAAlongMin+SAAshift1)) and (satlong < (SAAlongMax+SAAshift2))):
       SAAflag=1
       if (newSAA < 1):
          newSAA=1
          SAAbegin=elapsed
    else:
       if (newSAA == 1):
          SAAend=elapsed
          newSAA=0
       

    xm1=xmoon-xsat
    ym1=ymoon-ysat
    zm1=zmoon-zsat
    rm1=sqrt(xm1*xm1+ym1*ym1+zm1*zm1)

    moonangle=(acos((xm1*xsrc+ym1*ysrc+zm1*zsrc)/rm1)/dr)%360
    ramangle=(acos((vx*xsrc+vy*ysrc+vz*zsrc)/satspeed)/dr)%360
    earthangle=acos((-xsrc*xsat-ysrc*ysat-zsrc*zsat)/satdist)/dr
    Sun2Earth=acos((-xsun*xsat-ysun*ysat-zsun*zsat)/(satdist*rsun))/dr
    EarthConeAngle=asin(Rearth/satdist)/dr
    elevation=earthangle-EarthConeAngle
    SunElev=Sun2Earth-EarthConeAngle
# The Sun is considered a disc of 0.25 deg radius
    if SunElev < -0.25 :
       Nightflag=1
       if (newNight < 1):
          newNight=1
          Nightbegin=elapsed
    else:
       Nightflag=0
       if (newNight == 1):
          newNight=0
          Nightend=elapsed


    visible=0
    if (ramangle<minram):
       minram=ramangle
    if ((sunangle>SunAvoid) and (elevation>LimbAvoid)):
       if ((sunangle<(180-NegRollSunAvoid))):
          if ((moonangle>MoonAvoid)):
             if (SAAflag==0):
                lxczvis=lxczvis+SecInterval
                if (Nightflag==1):
                   sxtvis=sxtvis+SecInterval
                   uvtvis=uvtvis+SecInterval
                if (ramangle>RamAvoid):
                   visible=1
       
    Systime=(obsday-3654.0)*86400.0

    if (PrintOrb == 1):
       outfile.write('%12.2f %10.3f %10.3f %10.3f %7.3f %7.3f %7.3f %7.2f %6.2f %7.2f %7.2f %7.2f %7.2f %1d %1d %1d\n' % (Systime,xsat,ysat,zsat,vx,vy,vz,satlong,satlat,sunangle,moonangle,ramangle,elevation,SAAflag,Nightflag,visible))

  if (PrintOrb == 1):
     outfile.close()

  if (PrintVis == 1):
     visfile.close()

infile.close()

print ("")
