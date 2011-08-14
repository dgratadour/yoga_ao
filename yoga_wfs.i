require,"yoga_aolib.i";
require,"yoga_ao_ystruct.i";
require,"yoga_ao_utils.i"

func init_wfs_size(wfs,&pdiam,&Nfft,&Ntot,&nrebin,&pixsize,&qpixsize)
/* DOCUMENT init_wfs_size
   init_wfs_size,wfs,pdiam,Nfft,Ntot,nrebin,pixsize,qpixsize

   computes the quatum pixel sizes and all useful dimensions for a given wfs
   requires 2 externals : y_atmos & y_tel
   y_atmos : a atmos_struct()
   y_tel   : a tel_struct()
   wfs     : wfs_struct as input (input)
   pdiam   : pupil diam for each subap (pixels) (output)
   Nfft    : fft size for a subap (pixels) (output)
   Ntot    : hr image size for a subap (pixels) (output)
   nrebin  : rebin factor for a subap (output)
   pixsize : simulated pixel size for a subap (arcsec) (output)
   qpixsize: quantum pixel size for a subap (arcsec) (output)
   
   SEE ALSO:
 */
{
  /*
    Scheme to determine arrays sizes ...
    k = 5
    p = k * d/r0
    n = long(2*d*v/lambda/RASC)+1
    N = fft_goodsize(k*n/v*lambda/r0*RASC)
    u = k * lambda / r0 * RASC / N
    n = v/u - long(v/u) > 0.5 ? long(v/u)+1 : long(v/u)
    v = n * u
    Nt = v * Npix
   */

  extern y_atmos,y_tel;

  r0 = y_atmos.r0 / ((wfs.lambda/0.5)^(6./5));
  write,format="r0 for WFS : %.2f m\n",r0;
  seeing = RASC * (wfs.lambda * 1.e-6) / r0;   
  write,format="seeing for WFS : %.2f \"\n",seeing;

  if (pdiam == []) {
    // this case is usualy for the wfs with max # of subaps
    // we look for the best compromise between pixsize and fov
    
    subapdiam = y_tel.diam / double(wfs.nxsub); // diam of subap
    
    k = 6;
    pdiam = long(k * subapdiam / r0); // number of phase points per subap
    
    nrebin = long(2 * subapdiam * wfs.pixsize / (wfs.lambda*1.e-6) / RASC) + 1;
    // first atempt on a rebin factor
    Nfft = fft_goodsize(k * nrebin / wfs.pixsize * (wfs.lambda*1.e-6) / r0 * RASC);
    // size of the support in fourier domain
    
    qpixsize = k * (wfs.lambda*1.e-6) / r0 * RASC / Nfft;
    // quantum pixel size
  } else {
    // this case is for a wfs with fixed # of phase points
    subapdiam = y_tel.diam / double(wfs.nxsub); // diam of subap

    Nfft = fft_goodsize(2*pdiam);
    // size of the support in fourier domain

    qpixsize = pdiam * (wfs.lambda*1.e-6) / subapdiam * RASC / Nfft;
    // quantum pixel size
  }
  
  // actual rebin factor
  nrebin = ((wfs.pixsize/qpixsize - long(wfs.pixsize/qpixsize) > 0.5)(1) ?
            long(wfs.pixsize/qpixsize)+1 : long(wfs.pixsize/qpixsize));
  
  // actual pixel size
  pixsize = nrebin * qpixsize;

  if ((pixsize * wfs.npix > qpixsize * Nfft)(1)) 
    Ntot = long(pixsize * wfs.npix / qpixsize) + 1;
  else Ntot = Nfft;

  write,format="quantum pixsize : %.4f \"\n",qpixsize;
  write,format="simulated FoV : %.2f\" x %.2f\"\n",Ntot * qpixsize,Ntot * qpixsize;
  write,format="actual pixsize : %.2f\"\n",pixsize;
  write,format="actual FoV : %.2f\" x %.2f\"\n",pixsize * wfs.npix,
    pixsize * wfs.npix;
  write,format="number of phase points : %d\n",pdiam;
  write,format="size of fft support : %d\n",Nfft;
  
}

func init_wfs_geom(n,pupil)
/* DOCUMENT init_wfs_geom
   init_wfs_geom,n,pupil

   inits a wfs geometry (arrays used for image manipulation in the wfs model)
   requires 2 externals : y_atmos & y_wfs
   y_atmos : a atmos_struct()
   y_wfs   : a wfs_struct()
   n       : index of wfs to init
   pupil   : the pupil array
   
   SEE ALSO:
 */
{
  extern y_wfs;

  write,format="\n*-----------------------\nDoing inits on WFS # %d\n",n;

  if (n != wheremax(y_wfs.nxsub)(1))
    pdiam = (y_geom.pupdiam % y_wfs(n).nxsub == 0 ? long(y_geom.pupdiam / y_wfs(n).nxsub) :
             long(y_geom.pupdiam / y_wfs(n).nxsub) + 1);
  else pdiam = [];
  
  init_wfs_size,y_wfs(n),pdiam,Nfft,Ntot,nrebin,pixsize,qpixsize;
  
  y_wfs(n).pixsize   = pixsize;
  y_wfs(n)._pdiam    = pdiam;
  y_wfs(n)._Nfft     = Nfft;
  y_wfs(n)._Ntot     = Ntot;
  y_wfs(n)._nrebin   = nrebin;
  y_wfs(n)._qpixsize = qpixsize;
  y_wfs(n)._subapd   = y_tel.diam/y_wfs(n).nxsub;
  
  if (n == wheremax(y_wfs.nxsub)(1)) {
    //this is the wfs with largest # of subaps
    //the overall geometry is deduced from it
    geom_init,pdiam * y_wfs(n).nxsub;
  }
  
  // this is the i,j index of lower left pixel of subap
  istart = jstart = long(span(0.5,y_geom.pupdiam + 0.5,y_wfs(n).nxsub+1)+1)(:-1);
  y_wfs(n)._istart   = &istart;
  y_wfs(n)._jstart   = &jstart;

  // sorting out valid subaps
  fluxPerSub = array(float,y_wfs(n).nxsub,y_wfs(n).nxsub);

  for (i=1;i<=y_wfs(n).nxsub;i++) {
    for (j=1;j<=y_wfs(n).nxsub;j++) {
      indi = istart(i)+2;
      indj = jstart(j)+2;
      fluxPerSub(i,j) = sum((*y_geom._mpupil)(indi:indi+pdiam-1,indj:indj+pdiam-1));
    }    
  }
  fluxPerSub = fluxPerSub/pdiam^2.;
  tmp = fluxPerSub > y_wfs(n).fracsub;
  //tmp(where(tmp == 0)) = nvalid+10;
  y_wfs(n)._validsubs  = &int(where2(tmp)); 
  y_wfs(n)._isvalid    = &int(tmp);
  y_wfs(n)._nvalid     = sum(*y_wfs(n)._isvalid);
  y_wfs(n)._fluxPerSub = &float(fluxPerSub);

  
  // this defines how we cut the phase into subaps
  phasemap = array(0,pdiam,pdiam,y_wfs(n)._nvalid);
  tmp = indices(y_geom._n)-1; // we need c-like indices
  tmp = tmp(,,1)+tmp(,,2)*(y_geom._n);
  for (i=1;i<=y_wfs(n)._nvalid;i++) {
    indi = istart((*y_wfs(n)._validsubs)(1,i))+2;
    indj = jstart((*y_wfs(n)._validsubs)(2,i))+2;
    phasemap(,,i) = tmp(indi:indi+pdiam-1,indj:indj+pdiam-1);
  }
  //verif
  //pli,reform((*y_geom._mpupil)(*)(phasemap(*,1)+1),pdiam,pdiam)
  y_wfs(n)._phasemap = &int(phasemap(*,));

  //this is a phase shift of 1/2 pix in x and y
  halfxy = (1i * span(0,2*pi,y_wfs(n)._Nfft+1)(1:y_wfs(n)._pdiam) / 2.)(,-:1:y_wfs(n)._pdiam);
  halfxy += transpose(halfxy);

  //this defines how we create a larger fov if required
  if (y_wfs(n)._Ntot != y_wfs(n)._Nfft) {
    if (y_wfs(n)._Nfft % 2 != y_wfs(n)._Ntot % 2) 
      x1=long((y_wfs(n)._Ntot-y_wfs(n)._Nfft)/2.)+2;
    else
      x1=long((y_wfs(n)._Ntot-y_wfs(n)._Nfft)/2.)+1;
    x2=long(x1+y_wfs(n)._Nfft-1);

    tmp   = indices(y_wfs(n)._Nfft);
    hrpix = array(0.,y_wfs(n)._Ntot,y_wfs(n)._Ntot);
    hrpix(x1:x2,x1:x2) = roll(tmp(,,1)+(tmp(,,2)-1)*y_wfs(n)._Nfft);
    hrmap =where(roll(hrpix));
    
    y_wfs(n)._hrmap = &int(hrmap-1);
  } else y_wfs(n)._hrmap = &int(0);
  
  //creating the binindices array
  // note: if Ntot and nrebin*npix does not have the same parity
  // then it is impossible to rebin correctly while keeping
  // the center where it belongs: between 4 pix when Ntot is even
  // and on a pix when Ntot is odd
  // the +0.5 implies that when they don't have the same parity
  // the binned image is shifted by [-0.5,-0.5]
  // hence need to add half pixel tilt in x and y to the phase
  
  if (y_wfs(n)._nrebin*y_wfs(n).npix % 2 != y_wfs(n)._Ntot % 2) 
    x1=long((y_wfs(n)._Ntot-y_wfs(n)._nrebin*y_wfs(n).npix)/2.)+2;
  else
    x1=long((y_wfs(n)._Ntot-y_wfs(n)._nrebin*y_wfs(n).npix)/2.)+1;
  x2=long(x1+y_wfs(n)._nrebin*y_wfs(n).npix-1);
  
  if (y_wfs(n)._nrebin*y_wfs(n).npix % 2 != y_wfs(n)._Ntot % 2) 
    y_wfs(n)._halfxy = &float(halfxy.im);
  else y_wfs(n)._halfxy = &float(halfxy.im*0.);
  
  binindices = array(0,y_wfs(n)._Ntot,y_wfs(n)._Ntot);
  tmp = long((indices(y_wfs(n)._nrebin*y_wfs(n).npix) -1) / y_wfs(n)._nrebin);
  binindices(x1:x2,x1:x2) = tmp(,,1)+tmp(,,2)*y_wfs(n).npix+1;

  binmap = array(0,y_wfs(n)._nrebin*y_wfs(n)._nrebin,y_wfs(n).npix*y_wfs(n).npix);
  tmp = indices(y_wfs(n)._Ntot)-1;
  tmp = tmp(,,1)+tmp(,,2)*y_wfs(n)._Ntot;
  for (i=1;i<=y_wfs(n).npix*y_wfs(n).npix;i++) {
    if (y_wfs(n).gsalt > 0) binmap(,i) = tmp(where(binindices == i));
    else
      binmap(,i) = tmp(where(roll(binindices) == i));
  }
  y_wfs(n)._binmap = &int(binmap);

  /* verif
  fim = array(0.0f,y_wfs(n).npix,y_wfs(n).npix);
  for(i=1;i<=y_wfs(n).npix*y_wfs(n).npix;i++) fim(*)(i) = res(,,7)(*)(binmap(,i)+1)(sum);
  */
  
  //mapping the subaps binned images on the total image
  imamap =  array(0,y_wfs(n).npix*y_wfs(n).nxsub,y_wfs(n).npix*y_wfs(n).nxsub); 
  tmp = (indices(y_wfs(n).npix) -1);
  tmp = tmp(,,1)+tmp(,,2)*y_wfs(n).npix;
  for (i=1;i<=y_wfs(n)._nvalid;i++) {
    indi = ((*y_wfs(n)._validsubs)(1,i)-1)*(y_wfs(n).npix)+1;
    indj = ((*y_wfs(n)._validsubs)(2,i)-1)*(y_wfs(n).npix)+1;
    imamap(indi:indi+y_wfs(n).npix-1,indj:indj+y_wfs(n).npix-1) = tmp + (i-1) * (y_wfs(n).npix)^2 + 1;
  }
  y_wfs(n)._imamap = &int(where(imamap)(sort(imamap(where(imamap))))-1);
  
  // dealing with photometry
  telSurf  = pi/4.*y_tel.diam^2.;

  // from the guide star 
  if (y_wfs(n).gsalt == 0) {
    if (y_wfs(n).zerop == 0) y_wfs(n).zerop = 1e11;
    y_wfs(n)._nphotons = y_wfs(n).zerop*10^(-0.4*y_wfs(n).gsmag)*
      y_wfs(n).optthroughput*                 // include throughput to WFS
      (y_tel.diam/y_wfs(n).nxsub)^2./telSurf* // for unobstructed subaperture
      y_loop.ittime;                           // per iteration
  } else {  // we are dealing with a LGS
    y_wfs(n)._nphotons = y_wfs(n).lgsreturnperwatt*  // detected by WFS
      y_wfs(n).laserpower*                     // ... for given power
      y_wfs(n).optthroughput*                  // include throughput to WFS
      (y_tel.diam/y_wfs(n).nxsub)^2.*1e4*      // for unobstructed subaperture
      y_loop.ittime;                            // per iteration
  }
  write,format="nphotons : %.1f\n",y_wfs(n)._nphotons;
}


func make_lgs(proftype)
{
  extern y_wfs;

  //prof =  exp(-(altitude-9.e4)^2/(2*(4.e3)^2))*135.; // mono-modal
  // exp(-(altitude-8.7e4)^2/(2*(2.e3)^2))*110+exp(-(altitude-9.2e4)^2/(2*(2.e3)^2))*130; // bi-modal
  // exp(-(altitude-8.5e4)^2/(2*(1.5e3)^2))*55+exp(-(altitude-8.8e4)^2/(2*(2.e3)^2))*80+exp(-(altitude-9.2e4)^2/(2*(2.e3)^2))*120; // tri-modal
  // exp(-(altitude-8.7e4)^2/(2*(2.e3)^2))*130+exp(-(altitude-9.2e4)^2/(2*(2.e3)^2))*130; //bi-modal sym
  // note : the fwhm is 2*sqrt(2*log(2))*sigma ...


}

func prep_lgs_prof(numwfs,prof,h,beam,center=)
/* DOCUMENT prep_lgs_prof
   prep_lgs_prof,numwfs,prof,h,beam,center=

   The function returns an image array(double,n,n) of a laser beacon elongated by perpective
   effect. It is obtaind by convolution of a gaussian of width "lgsWidth" arcseconds, with the
   line of the sodium profile "prof". The altitude of the profile is the array "h".
   n        : image size
   z        : pixel size of the image (arcsec)
   lgsWidth : width of gaussian, in arcsec
   prof     : Na profile intensity, in arbitrary units
   h        : altitude, in meters. h MUST be an array with EQUALLY spaced elements.
   dOffAxis : offaxis distance of the subaperture that sees the beacon
   H        : altitude at which the beacon is supposed to be (in meters)
   azimut   : rotation of the beacon
   center   : string, either "image" or "fourier" depending on where the centre should be.
   
   Computation of LGS spot from the sodium profile:
   Everything is done here in 1D, because the Na profile is the result of the convolution of a function
   P(x,y) = profile(x) . dirac(y)
   by a gaussian function, for which variables x and y can be split :
   exp(-(x^2+y^2)/2.s^2)  =  exp(-x^2/2.s^2) * exp(-y^2/2.s^2)
   The convolution is (symbol $ denotes integral)
   C(X,Y) = $$ exp(-x^2/2.s^2) * exp(-y^2/2.s^2) * profile(x-X) * dirac(y-Y)  dx  dy
   First one performs the integration along y
   C(X,Y) = exp(-Y^2/2.s^2)  $ exp(-x^2/2.s^2) * profile(x-X)  dx
   which shows that the profile can be computed by
   - convolving the 1-D profile
   - multiplying it in the 2nd dimension by a gaussian function
     
   SEE ALSO:
 */
{
  extern y_wfs;
  
  subapdiam = y_tel.diam / double(y_wfs(numwfs).nxsub); // diam of subap
  if (y_wfs(numwfs).nxsub > 1)
    xsubs = span(-y_tel.diam/2+subapdiam/2,y_tel.diam/2-subapdiam/2,y_wfs(numwfs).nxsub);
  else xsubs = 0.;
  ysubs = xsubs;
  
  np = dimsof(prof)(2);       // number of points of the profile
  hG = sum(h*prof)/sum(prof); // center of gravity of the profile
  x=(indgen(y_wfs(numwfs)._Ntot)-(y_wfs(numwfs)._Ntot/2+1));
  // x expressed in pixels. (0,0) is in the fourier-center.
  x = x*y_wfs(numwfs)._qpixsize;  // x expressed in arcseconds
  dx = x(2)-x(1);

  
  if (y_wfs(numwfs).nxsub > 1)
    dOffAxis = sqrt((xsubs((*y_wfs(numwfs)._validsubs)(1,))-y_wfs(numwfs).lltx)^2 +
                    (ysubs((*y_wfs(numwfs)._validsubs)(2,))-y_wfs(numwfs).llty)^2);
  else 
    dOffAxis = sqrt((xsubs-y_wfs(numwfs).lltx)^2 +
                    (ysubs-y_wfs(numwfs).llty)^2);
  
  subsdone = array(1,y_wfs(numwfs)._nvalid);
  
  profi = array(0.0f,y_wfs(numwfs)._Ntot,y_wfs(numwfs)._nvalid);

  while (subsdone(*)(sum) > 0) {
    tmp = dOffAxis(where(subsdone)(1));
    inds = where(dOffAxis == tmp);
    // height, translated in arcsec due to perspective effect
    zhc = (h-hG)*(206265.*tmp/double(hG)^2.); 
    dzhc = zhc(2)-zhc(1);
    
    if (y_wfs(numwfs)._qpixsize > dzhc) {
      x_tmp = x(pcen);
      zhc_tmp = zhc(pcen);
      prof_tmp = prof(cum);
      profi(,inds) = interp(prof_tmp,zhc_tmp,x_tmp)(dif)*dzhc/dx;
    } else {
      x_tmp = x;
      zhc_tmp = zhc;
      prof_tmp = prof;
      profi(,inds) = interp(prof_tmp,zhc_tmp,x_tmp)*dzhc/dx;
    }
    subsdone(inds) = 0;
  }  
  
  /*
    place for manual interp
   */
  w = beam / 2.35482005;      //  sigma.   (2*sqrt(2*log(2)))=2.35482005
  if( w==0 ) {
    g = array(0.0,n);
    if( center=="image" )
      g(n/2:n/2+1)=0.5;
    else
      g(n/2+1)=1;
  }
  else {
    if( center=="image" )
      g = exp( -((x+z/2)^2/(2*w^2.) ) );
    else
      g = exp( -(x^2/(2*w^2.) ) );
    //=mygauss2(24,13,13,wx*2*sqrt(2*log(2)),wy*2*sqrt(2*log(2)),1.,0.,0.,,deriv=0);
  }

  // convolved profile in 1D.
  p1d = roll(fft(fft(profi,[1,0]) * fft(g(,-:1:y_wfs(numwfs)._nvalid),[1,0]),[-1,0]).re,[y_wfs(numwfs)._Ntot/2,0]);
  // abs is just to ensure only positive values (else values ~ -1e-12 may appear)
  p1d = abs(p1d);
  im = p1d(,-,) * g(-,,-:1:y_wfs(numwfs)._nvalid);
  if (center == "image")
    xcent = ycent = y_wfs(numwfs)._Ntot/2+0.5;
  else
    xcent = ycent = y_wfs(numwfs)._Ntot/2+1;

  for (cc=1;cc<=y_wfs(numwfs)._nvalid;cc++) {
    xsub = xsubs((*y_wfs(numwfs)._validsubs)(1,cc));
    ysub = ysubs((*y_wfs(numwfs)._validsubs)(2,cc));
    azimut = atan(ysub-y_wfs(numwfs).llty,xsub-y_wfs(numwfs).lltx);
    azimut = azimut*180./pi;
    im(,,cc) = fft_rotate(im(,,cc),-azimut);
    im(,,cc) /= sum(im(,,cc));
  }

  y_wfs(numwfs)._lgskern = &float(im);
} 


func wfs_map(arr,wfs,type=)
/* DOCUMENT wfs_map
   wfs_map,arr,wfs,type=

   maps an array of images onto a wfs
   arr     : the array to map
   wfs     : the wfs on which to map
   type    : type of mapping
   
   SEE ALSO:
 */
{
  if (type == []) type = "subaps"
  if (type == "subaps") {
    if (numberof(arr) != (*wfs._isvalid)(*)(sum))
      error,"wfs_map : wrong input dims";
    tmp = array(structof(arr),wfs.nxsub,wfs.nxsub);
    tmp(where(*wfs._isvalid)) = arr;
    return tmp;
  }
  if (type == "image") {
    if (dimsof(arr)(1) != 3)
      error,"wfs_map : wrong input dims";
    sz = dimsof(arr)(2);
    tmp = array(structof(arr),wfs.nxsub*sz,wfs.nxsub*sz);
    for (cc=1;cc<=wfs._nvalid;cc++) {
      indi = ((*wfs._validsubs)(1,cc)-1)*sz+1;
      indj = ((*wfs._validsubs)(2,cc)-1)*sz+1;
      tmp(indi:indi+sz-1,indj:indj+sz-1) = arr(,,cc);
    }
    return tmp;
  }
}

/*  //manual interp
  xprofbin = deltaprofbin = array(0,y_wfs(numwfs)._Ntot+1,y_wfs(numwfs)._nvalid);
    for (i=1;i<=numberof(x_tmp);i++) {
      if (x_tmp(i) > max(zhc_tmp)) xprofbin(i,cc) = numberof(zhc_tmp);
      else {
        tmp = where(zhc_tmp < x_tmp(i));
        if (numberof(tmp) != 0) {
          xprofbin(i,cc) = tmp(0);
          deltaprofbin(i,cc) = x_tmp(i) - zhc_tmp(tmp(0));
        } else xprofbin(i,cc) = 1;
      }
    }
    

    //interpolated :
    if (y_wfs(numwfs)._qpixsize > dzhc) {
      tmp = (prof_tmp(xprofbin(,cc)) + deltaprofbin(,cc)*(prof_tmp(clip(xprofbin(,cc)+1,,numberof(zhc_tmp)))-prof_tmp(xprofbin(,cc)))/dzhc);
      profi(,cc) = tmp(dif)*dzhc/dx;
    } else {
      tmp = (prof_tmp(xprofbin(1:y_wfs(numwfs)._Ntot,cc)) + deltaprofbin(1:y_wfs(numwfs)._Ntot,cc)*(prof_tmp(clip(xprofbin(1:y_wfs(numwfs)._Ntot,cc)+1,,numberof(zhc_tmp)))-prof_tmp(xprofbin(1:y_wfs(numwfs)._Ntot,cc)))/dzhc);
      profi(,cc) = tmp*dzhc/dx;
    }

 */
